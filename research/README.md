# Transit robustness experiments (2026-11-04)

This standalone research module is **synthetic**, separate from the published BLS demo and figures. It implements (1) stationary Gaussian AR(1) noise generation with a reproducible seed; (2) conditional box-depth GLS with AR(1) whitening, estimating baseline and depth jointly; and (3) a cadence-window geometric probability of observing at least three transit events over uniformly distributed epochs.

The GLS assumes equally spaced samples and known marginal noise and correlation. **Do not apply it directly across observing gaps**: split independent segments or construct a gap-aware covariance. The window statistic is not detection completeness: it ignores noise, BLS search efficiency and astrophysical false positives. No observational data or existing figures are changed.

Reproduce: `PYTHONPATH=research python -m pytest -q research/test_robustness.py`. Requires numpy and pytest. The tests check stationary moments, deterministic noise, exact noiseless depth recovery, correlation sensitivity, and missing-transit windows. This module is an experimental addition pending scientific review.
