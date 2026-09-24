"""Layer: data
WP: WP1, WP2, WP8, WP12
Claims: -
Status: implemented
Requirements: -

Content-addressed download cache and a polite HTTP client.

Every external file the simulator consumes lands under
``data/external/<source>/`` with a ``<file>.meta.json`` sidecar beside it.  The
sidecar records the URL it came from, the UTC time of retrieval, the server's
ETag or Last-Modified if either was offered, the SHA-256 of the bytes, and the
byte count.  That is the minimum needed for
:class:`tss.provenance.Manifest` to state which version of an upstream data set
produced a figure, and for a reader to tell whether their copy is the same one.

The filename is content-addressed in the weak sense that it carries a hash of
the request, not of the response.  Hashing the response would mean a file has
no name until it is downloaded, which makes resume and cache lookup awkward.
The response hash is recorded in the sidecar instead, and
:func:`verify_cached` re-checks it.

Politeness is not decoration.  These are public archives run on other people's
budgets, and one of them (wspr.live) explicitly asks that bulk users take the
published archive files rather than hammer the live endpoint.
:func:`polite_get` therefore sets a User-Agent that names the project and the
institution, enforces a minimum interval between requests to the same host,
retries 5xx with exponential backoff, and does not retry 4xx at all.  A 404 is
an answer, and asking again will not change it.

TLS verification is never disabled and proxy environment variables are never
unset.  If a certificate fails, that is a fact about the network which needs
fixing outside this package.
"""

from __future__ import annotations

import hashlib
import json
import os
import random
import re
import threading
import time
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import urlsplit

from tss.data.base import FetchError

__all__ = [
    "cache_root",
    "source_dir",
    "cached_path",
    "meta_path",
    "is_fresh",
    "record_meta",
    "read_meta",
    "verify_cached",
    "user_agent",
    "polite_get",
    "download",
    "DEFAULT_MIN_INTERVAL_S",
    "DEFAULT_TIMEOUT_S",
    "DEFAULT_MAX_RETRIES",
]


# ---------------------------------------------------------------------------
# Policy constants
# ---------------------------------------------------------------------------

DEFAULT_MIN_INTERVAL_S: float = 1.0
"""Minimum gap between two requests to the same host, seconds.

One second is not generous, it is merely not rude.  Individual modules raise it
where the provider asks for more.
"""

DEFAULT_TIMEOUT_S: float = 60.0
DEFAULT_MAX_RETRIES: int = 4
DEFAULT_BACKOFF_S: float = 2.0
DEFAULT_MAX_AGE_DAYS: float = 30.0

_RETRY_STATUS = frozenset({500, 502, 503, 504, 507, 509})
"""5xx codes worth retrying.  A 501 Not Implemented is not on this list, because
the server has told us the request will never work.
"""

_INSTITUTION = "Utah State University"
_PROJECT = "tss-sim (low-SWaP FMCW topside ionospheric sounder simulator)"

_last_request: dict[str, float] = {}
_throttle_lock = threading.Lock()

_SAFE = re.compile(r"[^A-Za-z0-9._-]+")


# ---------------------------------------------------------------------------
# Cache layout
# ---------------------------------------------------------------------------


def cache_root() -> Path:
    """The root of the external data cache.

    ``data/external/`` under the repository root, unless ``TSS_DATA_CACHE``
    overrides it.  Locating it from :func:`tss.config.repo_root` rather than
    from the current directory means a fetch started from anywhere puts its
    files in the same place, which is what lets the manifest paths mean
    something.
    """
    override = os.environ.get("TSS_DATA_CACHE")
    if override:
        return Path(override).expanduser()
    from tss.config import repo_root

    return repo_root() / "data" / "external"


def source_dir(source: str, *, create: bool = True) -> Path:
    """The cache directory for one source, created on demand."""
    d = cache_root() / _slug(source)
    if create:
        d.mkdir(parents=True, exist_ok=True)
    return d


def _slug(text: str, *, maxlen: int = 64) -> str:
    """A filesystem-safe fragment of ``text``, never empty."""
    out = _SAFE.sub("-", str(text)).strip("-.")
    if not out:
        out = "item"
    return out[:maxlen]


def cached_path(source: str, key: str, suffix: str = "") -> Path:
    """Where a cached item lives.

    ``key`` is whatever identifies the request: usually the URL, sometimes a
    URL plus a query string.  The returned name is a readable slug of the key's
    tail plus a short digest of the whole key, so that two requests that differ
    only in their query parameters do not collide and a human can still tell
    which file is which by looking at the directory.

    ``suffix`` is appended verbatim if the slug does not already end in it, so
    callers can ask for ``.asc``, ``.csv`` or ``.json`` and get a file whose
    extension tells the truth.
    """
    digest = hashlib.sha256(key.encode("utf-8")).hexdigest()[:12]
    tail = key.rstrip("/").rsplit("/", 1)[-1] or source
    tail = tail.split("?", 1)[0]
    stem = _slug(tail, maxlen=48) or source
    if suffix and stem.endswith(suffix):
        stem = stem[: -len(suffix)]
    return source_dir(source) / f"{stem}-{digest}{suffix}"


def meta_path(path: str | Path) -> Path:
    """The sidecar path for a cached file."""
    p = Path(path)
    return p.with_name(p.name + ".meta.json")


def is_fresh(path: str | Path, max_age_days: float = DEFAULT_MAX_AGE_DAYS) -> bool:
    """True if ``path`` exists, is non-empty, and is younger than the limit.

    A ``max_age_days`` that is None or not finite means "never stale", which is
    the right answer for a closed archive such as the ISIS topside files.  Age
    is taken from the sidecar's recorded retrieval time where one exists, and
    from the file mtime otherwise, because copying a tree preserves neither
    reliably but the sidecar is at least explicit.
    """
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        return False
    if max_age_days is None:
        return True
    try:
        if max_age_days != max_age_days or max_age_days == float("inf"):  # NaN or inf
            return True
    except TypeError:
        return True
    if max_age_days <= 0:
        return False

    retrieved = None
    meta = read_meta(p)
    if meta:
        stamp = meta.get("retrieved_utc")
        if isinstance(stamp, str):
            try:
                retrieved = datetime.fromisoformat(stamp)
            except ValueError:
                retrieved = None
    if retrieved is None:
        retrieved = datetime.fromtimestamp(p.stat().st_mtime, tz=timezone.utc)
    if retrieved.tzinfo is None:
        retrieved = retrieved.replace(tzinfo=timezone.utc)

    age_days = (datetime.now(timezone.utc) - retrieved).total_seconds() / 86400.0
    return age_days <= float(max_age_days)


# ---------------------------------------------------------------------------
# Sidecars
# ---------------------------------------------------------------------------


def record_meta(
    path: str | Path,
    *,
    url: str = "",
    etag: str | None = None,
    last_modified: str | None = None,
    source: str = "",
    extra: Mapping[str, Any] | None = None,
) -> Path:
    """Write the ``<file>.meta.json`` sidecar and return its path.

    The SHA-256 comes from :func:`tss.provenance.file_digest`, so the digest in
    a sidecar and the digest in a run manifest are computed by the same code
    and can be compared without qualification.
    """
    from tss.provenance import file_digest

    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"cannot record metadata for missing file {p}")

    payload: dict[str, Any] = {
        "url": url,
        "source": source,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "etag": etag,
        "last_modified": last_modified,
        "sha256": file_digest(p),
        "bytes": p.stat().st_size,
        "filename": p.name,
        "user_agent": user_agent(),
    }
    if extra:
        payload["extra"] = dict(extra)

    out = meta_path(p)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return out


def read_meta(path: str | Path) -> dict[str, Any] | None:
    """Read a sidecar, or None if there is none or it is unreadable.

    A corrupt sidecar returns None rather than raising.  The sidecar is a record
    about a file, not the file, and losing it should degrade provenance rather
    than stop a run.  Callers that need provenance to be present check for None.
    """
    mp = meta_path(path)
    if not mp.is_file():
        return None
    try:
        data = json.loads(mp.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return data if isinstance(data, dict) else None


def verify_cached(path: str | Path) -> bool:
    """True if the file's bytes still match the digest in its sidecar."""
    from tss.provenance import file_digest

    meta = read_meta(path)
    if not meta or not meta.get("sha256"):
        return False
    return file_digest(path) == meta["sha256"]


# ---------------------------------------------------------------------------
# HTTP
# ---------------------------------------------------------------------------


def user_agent() -> str:
    """The User-Agent string sent with every request.

    It names the project and the institution so that an archive operator who
    sees the traffic can tell who it is and why.  ``TSS_USER_AGENT`` overrides
    it for anyone running this from somewhere else.  No personal contact detail
    is placed in the header.
    """
    override = os.environ.get("TSS_USER_AGENT")
    if override:
        return override
    try:
        from tss import __version__ as version
    except Exception:  # pragma: no cover - version is cosmetic here
        version = "0"
    return f"tss-sim/{version} ({_PROJECT}; research use; {_INSTITUTION})"


def _throttle(url: str, min_interval_s: float) -> None:
    """Sleep until at least ``min_interval_s`` has passed for this host."""
    if min_interval_s <= 0:
        return
    host = urlsplit(url).netloc.lower()
    with _throttle_lock:
        previous = _last_request.get(host)
        now = time.monotonic()
        if previous is not None:
            wait = min_interval_s - (now - previous)
            if wait > 0:
                time.sleep(wait)
                now = time.monotonic()
        _last_request[host] = now


def polite_get(
    url: str,
    *,
    params: Mapping[str, Any] | None = None,
    headers: Mapping[str, str] | None = None,
    timeout: float = DEFAULT_TIMEOUT_S,
    max_retries: int = DEFAULT_MAX_RETRIES,
    min_interval_s: float = DEFAULT_MIN_INTERVAL_S,
    backoff_s: float = DEFAULT_BACKOFF_S,
    stream: bool = False,
    session: Any | None = None,
) -> Any:
    """GET a URL politely and return the :class:`requests.Response`.

    Retries are attempted only for the codes in ``_RETRY_STATUS`` and for
    transport-level errors.  A 4xx is returned to the caller as a
    :class:`~tss.data.base.FetchError` immediately and without a second
    request: a 400 means the query was wrong, a 403 means we are not allowed,
    and a 429 means we have already asked too often, so in every case the right
    move is to stop rather than to try again.  429 in particular is excluded
    from the retry set on purpose; honouring it means backing off for longer
    than a fetch should block for, and the caller should reschedule.

    ``params`` may be a mapping or a sequence of ``(key, value)`` pairs.  The
    pair form is what a service that takes a repeated key needs; Madrigal's
    ``isprintService.py`` takes one ``parms=`` per parameter, for instance.

    ``requests`` is imported here rather than at module scope so that importing
    :mod:`tss.data` costs nothing on a machine that will only read cached files.
    """
    import requests  # imported lazily, see the docstring

    sess = session if session is not None else requests

    def _send(hdrs: dict[str, str]) -> Any:
        # verify is left at its default.  Never pass verify=False here.
        return sess.get(
            url,
            params=_as_params(params),
            headers=hdrs,
            timeout=timeout,
            stream=stream,
        )

    return _polite_request(
        "GET", url, _send, headers=headers, max_retries=max_retries,
        min_interval_s=min_interval_s, backoff_s=backoff_s,
    )


def _as_params(params: Any) -> Any:
    """Pass a mapping through as a dict and a sequence of pairs as a list."""
    if not params:
        return None
    if isinstance(params, Mapping):
        return dict(params)
    return list(params)


def _polite_request(
    method: str,
    url: str,
    send: Any,
    *,
    headers: Mapping[str, str] | None,
    max_retries: int,
    min_interval_s: float,
    backoff_s: float,
) -> Any:
    """The retry loop behind :func:`polite_get`, kept separate so a second method
    could share it without duplicating the policy."""
    hdrs: dict[str, str] = {"User-Agent": user_agent(), "Accept-Encoding": "gzip"}
    if headers:
        hdrs.update(headers)

    last_error: Exception | None = None
    for attempt in range(max_retries + 1):
        _throttle(url, min_interval_s)
        try:
            response = send(hdrs)
        except Exception as exc:  # transport failure, worth one more try
            last_error = exc
            if attempt >= max_retries:
                raise FetchError(
                    f"{method} {url} failed after {attempt + 1} attempts: {exc}"
                ) from exc
            _sleep_backoff(backoff_s, attempt)
            continue

        status = int(getattr(response, "status_code", 0))
        if 200 <= status < 400:
            return response

        if 400 <= status < 500:
            reason = getattr(response, "reason", "") or ""
            body = _short_body(response)
            raise FetchError(
                f"{method} {url} returned {status} {reason}; not retrying a client error. "
                f"{body}".strip()
            )

        if status in _RETRY_STATUS and attempt < max_retries:
            _sleep_backoff(backoff_s, attempt)
            continue

        raise FetchError(
            f"{method} {url} returned {status}; giving up after {attempt + 1} attempts"
        )

    raise FetchError(f"{method} {url} failed: {last_error}")  # pragma: no cover


def _sleep_backoff(backoff_s: float, attempt: int) -> None:
    """Exponential backoff with a little jitter, so parallel runs desynchronise."""
    delay = backoff_s * (2.0**attempt)
    time.sleep(delay * (0.75 + 0.5 * random.random()))


def _short_body(response: Any, limit: int = 200) -> str:
    """A truncated response body, for putting the server's complaint in the error."""
    try:
        text = response.text
    except Exception:
        return ""
    if not text:
        return ""
    text = " ".join(str(text).split())
    return f"Server said: {text[:limit]}"


def download(
    url: str,
    dest: str | Path,
    *,
    source: str = "",
    force: bool = False,
    max_age_days: float = DEFAULT_MAX_AGE_DAYS,
    chunk: int = 1 << 16,
    extra_meta: Mapping[str, Any] | None = None,
    **get_kwargs: Any,
) -> Path:
    """Fetch ``url`` to ``dest``, writing a sidecar, and return ``dest``.

    A cached file that is fresh and whose digest still matches its sidecar is
    returned without a request.  The download goes to a ``.part`` file and is
    renamed on success, so an interrupted fetch never leaves a truncated file
    that looks complete.
    """
    dest = Path(dest)
    if not force and is_fresh(dest, max_age_days) and verify_cached(dest):
        return dest

    dest.parent.mkdir(parents=True, exist_ok=True)
    response = polite_get(url, stream=True, **get_kwargs)
    tmp = dest.with_name(dest.name + ".part")
    try:
        with open(tmp, "wb") as fh:
            for block in response.iter_content(chunk_size=chunk):
                if block:
                    fh.write(block)
        tmp.replace(dest)
    finally:
        if tmp.exists():
            tmp.unlink()

    headers = getattr(response, "headers", {}) or {}
    record_meta(
        dest,
        url=url,
        etag=headers.get("ETag"),
        last_modified=headers.get("Last-Modified"),
        source=source,
        extra=extra_meta,
    )
    return dest
