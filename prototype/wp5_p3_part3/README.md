# WP5 P3 part 3 prototype: NOT FOR MERGE

Reference only.  Written in a cloud session on 2026-11-27 to check approaches
before production code is written on Andrew's machine.  Andrew's rule
(2026-11-27): production code is written locally, so nothing on this branch
is merged.  The local session reads it and writes the code and tests fresh.

## What is here

1. `src/tss/l3_signal/drift.py`: the detector's false-alarm rate stated per
   sequence (decisions (m), (q)), and `estimator="joint"` (decisions (o), (r),
   (s), (t)).  Tests in `tests/unit/test_l3_drift.py` (the per-sequence tests
   passed in the cloud).
2. `src/tss/l3_signal/novel/route_one_joint.py`: the deterministic joint fit
   (variable projection, complex amplitude solved in closed form).
3. `prototype/wp5_p3_part3/fading_ml3.py`: the fading-aware maximum-likelihood
   fit, as an experiment script only (not yet in `src/`).

## What the experiments found (parameter level, the part 2 scene)

Deterministic, 200 seeds (`joint_mc.txt`): the joint fit is 1.02, 1.02 and
1.11 times the bound at 20, 15 and 10 dB; the pinned estimator 1.51, 1.72,
2.29.  Its sigma matches the spread.

Rayleigh fading (`joint_fading.txt`): pinned and deterministic joint are 4 to
9 times the fading bound (about 20 to 30 m/s).

Fading-aware ML (`fading_ml3.py`, 20 seeds, formed facets only): 0.85 to 1.16
times the fading bound at tau_c 0.5 and 2 s, 20 and 10 dB, with tau_c known or
fitted.  Two things were needed:
- the noise level fixed at the measured value (the tier has it from the
  negative-beat half); left free, bad fits explained a fading spike as a
  narrow facet plus extra noise;
- a coarse grid over (t0, drift, width) before the local fit; started from the
  pinned estimate alone, the local fit stopped in those minima.
It costs about 3 to 6 s per echo on one core in the cloud (grid of 1690
likelihood evaluations at N = 209), which matters for the study's run time.
