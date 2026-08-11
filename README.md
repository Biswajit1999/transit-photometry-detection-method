# Transit Photometry — Exoplanet Detection Method

How the majority of known exoplanets have been found: watching a star's
brightness for the tiny, periodic dip caused by a planet passing in
front of it. This repo explains the physics and statistics of the
method and implements a real Box Least Squares (BLS) period search —
the actual algorithm class the Kepler and TESS pipelines use — from
scratch in Python, validated by injecting a known signal and recovering
it.

## The physics

When a planet of radius $R_p$ transits a star of radius $R_\star$, it
blocks a fraction of the starlight equal to the ratio of their disk
areas:

$$\delta \approx \left(\frac{R_p}{R_\star}\right)^2$$

For an Earth-Sun-like pair this is about 84 parts per million; for a
Jupiter-Sun-like pair it is roughly 1%. The transit duration depends on
the orbital period $P$, the star's radius, the orbital semi-major axis
$a$, and the impact parameter $b$ (how centrally the planet crosses the
disk):

$$T_{dur} \approx \frac{P}{\pi} \arcsin\left(\frac{R_\star}{a}\sqrt{(1+R_p/R_\star)^2 - b^2}\right)$$

Three real transits at consistent period, depth, and duration are
usually required before a signal is called a genuine planet candidate,
since single dips can be caused by instrumental artifacts, stellar
activity, or eclipsing binaries in the background.

## Why this method dominates the real exoplanet census

Transit photometry needs only a stable, well-calibrated brightness
measurement — no need to resolve the star and planet separately, and no
need for the extreme spectral precision radial velocity requires. This
makes it well suited to wide-field surveys that monitor huge numbers of
stars simultaneously. As of 2026, the majority of the ~5,800+ confirmed
exoplanets in the NASA Exoplanet Archive were found this way, primarily
via the Kepler, K2, and TESS missions.

**Real limitation:** transit photometry only detects planets whose
orbits happen to be aligned edge-on as seen from Earth — a geometric
probability of roughly $R_\star/a$, meaning it systematically
undercounts real planetary systems whose orbits aren't so aligned, and
it directly measures $R_p$, not mass (radial velocity or transit-timing
follow-up is needed for that).

## What this repo's code does

`scripts/transit_photometry_demo.py`:

1. Injects a known transit signal (period, depth, duration) into a
   synthetic light curve with **Kepler's own published ~100 ppm per-30-
   minute photometric noise level** for a quiet Sun-like star — a real,
   published precision regime, not an arbitrary number.
2. Implements Box Least Squares (Kovacs, Zsom & Mazeh 2002) from
   scratch: a 2-D search over trial period *and* trial phase offset,
   sliding a box of the expected transit duration across the folded
   light curve and keeping the period/phase combination with the
   strongest signal residue — the same core logic used by the Kepler
   and TESS transit-search pipelines.
3. Recovers the period, depth, and detection significance, and reports
   the error against the known injected "ground truth" — a genuine
   validation of the method's real-world sensitivity.

Run it yourself:

```bash
pip install -r requirements.txt
python scripts/transit_photometry_demo.py
```

## Result

| Quantity | Injected | Recovered | Error |
|---|---|---|---|
| Period | 4.35 days | 4.3500 days | 0.00% |
| Depth | 850.0 ppm | 868.3 ppm | 2.15% |
| Detection significance | — | 74.2σ | — |

The BLS period search cleanly isolates the true period against strong
aliasing at harmonics/sub-harmonics (visible as smaller peaks in
`figures/transit_bls_recovery.png`), and the phase-folded light curve
shows a clean, textbook box-shaped transit recovered from noisy
individual measurements.

## Why this repo uses simulated (not raw archival) data

This repo demonstrates the *method itself* — its sensitivity, biases,
and failure modes — which is best shown with a known "ground truth" to
validate recovery against. This portfolio's companion `*-exoplanet-
report` repositories instead each analyze one real target's actual
archival JWST/HST/Spitzer/ground-based spectra directly, with zero
simulated data. Both approaches are stated plainly here rather than
blurring the two.

## Repository structure

```text
scripts/transit_photometry_demo.py   BLS implementation + injection-recovery test
figures/                             generated plot + summary_statistics.csv
```

## References

1. Charbonneau, D. et al., 2000. Detection of Planetary Transits Across
   a Sun-like Star. *The Astrophysical Journal Letters*, 529(1), L45 —
   the first exoplanet transit detection (HD 209458 b).
2. Kovacs, G., Zsom, A. and Mazeh, T., 2002. A box-fitting algorithm in
   the search for periodic transits. *Astronomy & Astrophysics*, 391,
   pp.369-377.
3. Mandel, K. and Agol, E., 2002. Analytic Light Curves for Planetary
   Transit Searches. *The Astrophysical Journal Letters*, 580(2), L171.
4. Borucki, W.J. et al., 2010. Kepler Planet-Detection Mission:
   Introduction and First Results. *Science*, 327(5968), pp.977-980.
5. Ricker, G.R. et al., 2015. Transiting Exoplanet Survey Satellite
   (TESS). *Journal of Astronomical Telescopes, Instruments, and
   Systems*, 1(1), 014003.
6. NASA Exoplanet Archive, <https://exoplanetarchive.ipac.caltech.edu/>.

## Author

Biswajit Jana — [Portfolio](https://biswajit1999.github.io/Biswajit_Jana.github.io/) · [GitHub](https://github.com/Biswajit1999) · [LinkedIn](https://www.linkedin.com/in/biswajit-jana-27011a151/) · [ORCID](https://orcid.org/0009-0002-2411-1891)
