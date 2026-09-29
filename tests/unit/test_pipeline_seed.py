"""The pipeline's random stream must not depend on Python's hash salt.

``Pipeline._rng`` once folded the configuration name in with the built-in
``hash``, which Python salts per process.  Every run of a campaign then drew
different noise from the same seed, and a scored result could not be
reproduced.  This test runs the draw in two child processes with different
``PYTHONHASHSEED`` values and requires the same numbers from both.
"""

from __future__ import annotations

import os
import subprocess
import sys

SNIPPET = (
    "from tss.pipeline import Pipeline\n"
    "print(repr(Pipeline({}, seed=7)._rng('leo_500km_single').random(4).tolist()))\n"
)


def _draw(hash_seed: str) -> str:
    env = dict(os.environ, PYTHONHASHSEED=hash_seed)
    done = subprocess.run([sys.executable, "-c", SNIPPET], env=env,
                          capture_output=True, text=True, check=True)
    return done.stdout.strip()


def test_rng_is_independent_of_the_hash_salt():
    assert _draw("1") == _draw("12345")


def test_rng_still_differs_between_configurations():
    from tss.pipeline import Pipeline

    p = Pipeline({}, seed=7)
    assert p._rng("a").random() != p._rng("b").random()
