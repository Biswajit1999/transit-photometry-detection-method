# Ephemeris and forecast sensitivity — 21 November 2026

Two conditional numerical diagnostics: weighted linear least-squares estimates of reference transit epoch and period, and propagation of their full 2×2 covariance into future timing uncertainty. Cycle numbers must be known integers (or explicitly defined epoch indices); midpoint and sigma must share time units. The covariance is formal under independent Gaussian timing errors and does not include correlated timing noise, period derivatives, TTVs, or ephemeris-source systematics.

Run from this folder: `python -m unittest discover -p 'test_*.py' -v`. Three synthetic tests were run against equivalent local development code and passed; the exact published GitHub branch has not been independently executed. This does not alter the existing BLS results or claim a new planet. AI-assisted code preparation disclosed; Biswajit Jana remains the requested sole research author. Review exact Git commit author metadata and CI before merging.
