"""Transit photometry method demonstration: inject a transit signal into
a synthetic light curve with Kepler-class photometric noise, then
recover the period, phase, duration, and depth using a Box Least
Squares (BLS) search (Kovacs, Zucker & Mazeh 2002) -- a standard, widely
used transit-search algorithm, related in purpose to but not the same
as the adaptive wavelet matched-filter architecture the actual
Kepler/TESS mission pipelines run. Duration is searched over a grid,
not handed to the algorithm, so the search is blind to all four
injected parameters, not just two of them.

This is a PEDAGOGICAL DEMONSTRATION with simulated data, not a specific
real target's raw archival light curve (see README.md for why, and see
this portfolio's *-exoplanet-report repos for 11 planets analyzed
directly from real archival JWST/HST/Spitzer/ground-based data). The
injected period/depth are in a hot-Neptune-class range broadly similar
to real Kepler discoveries, and the noise level matches Kepler's own
published ~100 ppm per-30-min photometric precision for a quiet
Sun-like star.

CAVEAT: the quoted "detection SNR" is the in-transit depth
signal-to-noise under the known Gaussian noise model at the recovered
period/phase/duration -- it is not a trial-corrected false-alarm
probability accounting for the number of period/phase/duration
combinations searched, which is what a real survey pipeline's
significance threshold has to account for. The synthetic light curve
also has none of the red noise, data gaps, stellar variability, or
limb-darkened ingress/egress shape a real Kepler or TESS light curve
would show, so this demo's recovery is easier than a real low-SNR
candidate search.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update(
    {
        "axes.spines.right": False,
        "axes.spines.top": False,
        "figure.facecolor": "white",
        "font.size": 9,
        "savefig.bbox": "tight",
    }
)

FIG_DIR = Path(__file__).resolve().parents[1] / "figures"

# Injected "ground truth" transit parameters (realistic hot-Neptune regime,
# broadly similar to real Kepler-class discoveries -- not a specific target).
TRUE_PERIOD_DAYS = 4.35
TRUE_DEPTH_PPM = 850.0
TRUE_DURATION_HR = 2.8
TRUE_T0_DAYS = 1.7

BASELINE_DAYS = 60.0
CADENCE_MIN = 30.0  # Kepler long-cadence
NOISE_PPM = 100.0  # real Kepler-class per-30-min photometric precision for a quiet star
N_PHASE_BINS = 200
SEARCH_PERIODS = np.arange(1.0, 15.0, 0.025)
NULL_TRIALS = 100
RECOVERY_TRIALS = 20
RECOVERY_DEPTHS_PPM = (25.0, 50.0, 75.0, 100.0, 150.0, 300.0, 850.0)

# Trial transit durations to search, as a fraction of the trial period --
# spans a wider range than the injected 2.8 h so the search doesn't just
# get handed the answer. Real BLS implementations (e.g. the NASA Exoplanet
# Archive's periodogram service) search a comparable fractional-duration
# range rather than a single fixed duration.
DURATION_FRACTIONS = np.array([0.01, 0.02, 0.03, 0.05, 0.08, 0.12])


def make_cadence(rng: np.random.Generator) -> np.ndarray:
    """Return a frozen long-cadence window with declared interruptions."""
    time = np.arange(0, BASELINE_DAYS, CADENCE_MIN / (24 * 60))
    scheduled = np.ones(time.size, dtype=bool)
    for gap_start in (14.5, 29.5, 44.5):
        scheduled &= ~((time >= gap_start) & (time < gap_start + 0.75))
    random_keep = rng.random(time.size) >= 0.03
    return time[scheduled & random_keep]


def transit_model(
    time: np.ndarray,
    period: float,
    t0: float,
    depth_ppm: float,
    duration_hr: float,
) -> np.ndarray:
    phase = ((time - t0 + period / 2) % period) - period / 2
    half_dur_days = (duration_hr / 24.0) / 2.0
    taper = half_dur_days * 0.15  # smoothed ingress/egress
    flux = np.ones_like(time)
    in_transit = np.abs(phase) < half_dur_days
    edge = (np.abs(phase) >= half_dur_days) & (np.abs(phase) < half_dur_days + taper)
    flux[in_transit] -= depth_ppm * 1e-6
    edge_frac = 1.0 - (np.abs(phase[edge]) - half_dur_days) / taper
    flux[edge] -= depth_ppm * 1e-6 * edge_frac
    return flux


def best_box_for_period(
    time: np.ndarray,
    flux: np.ndarray,
    period: float,
    duration_hr: float,
    n_phase_bins: int = N_PHASE_BINS,
) -> dict[str, float | int]:
    """Bin the folded light curve in phase and slide a box of the expected
    duration across ALL phase offsets (with wraparound), returning the
    offset and stats of the best-fitting box. Scanning every offset -- not
    just one fixed reference phase -- is essential: a real transit can fall
    anywhere in phase. This is the two-dimensional (period x phase) search
    that Kovacs, Zucker & Mazeh (2002)'s BLS performs. Kepler's
    production search used a distinct wavelet-based matched filter.
    """
    phase = (time % period) / period
    bin_idx = np.clip((phase * n_phase_bins).astype(int), 0, n_phase_bins - 1)
    bin_sum = np.bincount(bin_idx, weights=flux, minlength=n_phase_bins)
    bin_count = np.bincount(bin_idx, minlength=n_phase_bins)
    total_sum, total_count = bin_sum.sum(), bin_count.sum()

    box_width_bins = max(1, round((duration_hr / 24.0 / period) * n_phase_bins))
    ext_sum = np.concatenate([bin_sum, bin_sum[:box_width_bins]])
    ext_count = np.concatenate([bin_count, bin_count[:box_width_bins]])
    window_sum = np.convolve(ext_sum, np.ones(box_width_bins), mode="valid")[:n_phase_bins]
    window_count = np.convolve(ext_count, np.ones(box_width_bins), mode="valid")[:n_phase_bins]

    valid = window_count >= 3
    out_sum = total_sum - window_sum
    out_count = total_count - window_count
    with np.errstate(invalid="ignore", divide="ignore"):
        in_mean = window_sum / window_count
        out_mean = out_sum / np.maximum(out_count, 1)
        depth = out_mean - in_mean
        sr = depth * np.sqrt(window_count * out_count / np.maximum(window_count + out_count, 1))
    sr = np.where(valid, sr, -np.inf)

    best_bin = int(np.argmax(sr))
    return {
        "power": max(float(sr[best_bin]), 0.0),
        "phase_start": best_bin / n_phase_bins,
        "phase_width": box_width_bins / n_phase_bins,
        "depth": float(depth[best_bin]),
        "n_in": int(window_count[best_bin]),
    }


def best_box_multi_duration(
    time: np.ndarray,
    flux: np.ndarray,
    period: float,
    duration_fractions: np.ndarray,
) -> dict[str, float | int]:
    """Search a grid of trial durations (as a fraction of the trial
    period) rather than being told the true duration, and keep whichever
    gives the strongest signal residue -- this is what makes the search
    blind to the injected duration, not just the injected period/phase."""
    best = None
    for frac in duration_fractions:
        duration_hr = frac * period * 24.0
        result = best_box_for_period(time, flux, period, duration_hr)
        result["duration_hr"] = duration_hr
        if best is None or result["power"] > best["power"]:
            best = result
    return best


def box_least_squares(
    time: np.ndarray,
    flux: np.ndarray,
    periods: np.ndarray,
    duration_fractions: np.ndarray,
) -> np.ndarray:
    power = np.zeros_like(periods)
    for i, p in enumerate(periods):
        power[i] = best_box_multi_duration(time, flux, p, duration_fractions)["power"]
    return power


def wilson_interval(
    successes: int, trials: int, z: float = 1.959963984540054
) -> tuple[float, float]:
    fraction = successes / trials
    denominator = 1 + z**2 / trials
    centre = (fraction + z**2 / (2 * trials)) / denominator
    half_width = (
        z
        * np.sqrt(
            fraction * (1 - fraction) / trials + z**2 / (4 * trials**2)
        )
        / denominator
    )
    return max(0.0, centre - half_width), min(1.0, centre + half_width)


def search(time: np.ndarray, flux: np.ndarray) -> tuple[np.ndarray, float, dict]:
    """Run the declared period, phase, and duration search."""
    power = box_least_squares(time, flux, SEARCH_PERIODS, DURATION_FRACTIONS)
    best_index = int(np.argmax(power))
    best_period = float(SEARCH_PERIODS[best_index])
    best_box = best_box_multi_duration(
        time, flux, best_period, DURATION_FRACTIONS
    )
    return power, best_period, best_box


def calibrate_null(
    time: np.ndarray,
    rng: np.random.Generator,
    trials: int = NULL_TRIALS,
) -> np.ndarray:
    """Maximum BLS residue under white noise on the exact timestamps."""
    maxima = np.empty(trials)
    for index in range(trials):
        flux = 1 + rng.normal(0, NOISE_PPM * 1e-6, size=time.size)
        power, _, _ = search(time, flux)
        maxima[index] = power.max()
    return maxima


def recovery_experiment(
    time: np.ndarray,
    threshold: float,
    rng: np.random.Generator,
    trials: int = RECOVERY_TRIALS,
) -> list[dict[str, float | int | str]]:
    """Recover randomized-epoch injections on one fixed timestamp window."""
    rows = []
    for depth_ppm in RECOVERY_DEPTHS_PPM:
        successes = 0
        recovered_periods = []
        for _ in range(trials):
            epoch = rng.uniform(0, TRUE_PERIOD_DAYS)
            flux = transit_model(
                time,
                TRUE_PERIOD_DAYS,
                epoch,
                depth_ppm,
                TRUE_DURATION_HR,
            )
            flux += rng.normal(0, NOISE_PPM * 1e-6, size=time.size)
            power, best_period, _ = search(time, flux)
            recovered_periods.append(best_period)
            if (
                power.max() >= threshold
                and abs(best_period / TRUE_PERIOD_DAYS - 1) <= 0.01
            ):
                successes += 1
        low, high = wilson_interval(successes, trials)
        rows.append(
            {
                "depth_ppm": depth_ppm,
                "trials": trials,
                "recoveries": successes,
                "recovery_fraction": successes / trials,
                "wilson95_low": low,
                "wilson95_high": high,
                "median_best_period_days": float(np.median(recovered_periods)),
                "criterion": (
                    "max residue >= cadence null 99th percentile and "
                    "|Pbest/Pinj-1| <= 1%"
                ),
            }
        )
    return rows


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)
    rng = np.random.default_rng(seed=11)

    time = make_cadence(rng)
    flux = transit_model(time, TRUE_PERIOD_DAYS, TRUE_T0_DAYS, TRUE_DEPTH_PPM, TRUE_DURATION_HR)
    flux += rng.normal(0, NOISE_PPM * 1e-6, size=time.size)

    power, best_period, best_box = search(time, flux)
    recovered_duration_hr = best_box["duration_hr"]
    recovered_depth_ppm = best_box["depth"] * 1e6
    depth_err_ppm = NOISE_PPM / np.sqrt(best_box["n_in"])
    snr = recovered_depth_ppm / depth_err_ppm

    period_error_pct = abs(best_period - TRUE_PERIOD_DAYS) / TRUE_PERIOD_DAYS * 100
    duration_error_pct = abs(recovered_duration_hr - TRUE_DURATION_HR) / TRUE_DURATION_HR * 100
    depth_error_pct = abs(recovered_depth_ppm - TRUE_DEPTH_PPM) / TRUE_DEPTH_PPM * 100

    null_maxima = calibrate_null(time, rng)
    null_threshold_99 = float(np.quantile(null_maxima, 0.99, method="higher"))
    observed_peak = float(power.max())
    null_exceedances = int(np.sum(null_maxima >= observed_peak))
    fap_plus_one = (null_exceedances + 1) / (NULL_TRIALS + 1)
    recovery_rows = recovery_experiment(time, null_threshold_99, rng)

    # Fold using the recovered box center and recovered duration so the
    # diagnostic plot reflects what the search actually found, not the
    # injected ground truth.
    box_center_phase = best_box["phase_start"] + best_box["phase_width"] / 2
    phase_days = (((time / best_period - box_center_phase + 0.5) % 1.0) - 0.5) * best_period
    half_dur_days = (recovered_duration_hr / 24.0) / 2.0

    summary_path = FIG_DIR / "summary_statistics.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "injected", "recovered", "error_pct"])
        writer.writerow(["period_days", TRUE_PERIOD_DAYS, f"{best_period:.4f}", f"{period_error_pct:.2f}"])
        writer.writerow(["duration_hr", TRUE_DURATION_HR, f"{recovered_duration_hr:.2f}", f"{duration_error_pct:.2f}"])
        writer.writerow(["depth_ppm", TRUE_DEPTH_PPM, f"{recovered_depth_ppm:.1f}", f"{depth_error_pct:.2f}"])
        writer.writerow(["detection_snr_uncorrected", "-", f"{snr:.1f}", "-"])
        writer.writerow(["n_transits_in_baseline", int(BASELINE_DAYS // TRUE_PERIOD_DAYS), "-", "-"])
        writer.writerow(["cadence_points", int(BASELINE_DAYS * 48), len(time), "scheduled + random gaps"])
        writer.writerow(["cadence_duty_cycle", 1.0, f"{len(time) / (BASELINE_DAYS * 48):.5f}", "-"])
        writer.writerow(["bls_peak_residue", "noise-only null", f"{observed_peak:.10f}", "-"])
        writer.writerow(["null_99pct_max_residue", f"{NULL_TRIALS} white-noise trials", f"{null_threshold_99:.10f}", "-"])
        writer.writerow(["mc_false_alarm_exceedances", NULL_TRIALS, null_exceedances, "-"])
        writer.writerow(["mc_fap_plus_one", "cadence + known white noise", f"{fap_plus_one:.8f}", "-"])

    sensitivity_path = FIG_DIR / "depth_recovery.csv"
    with sensitivity_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=recovery_rows[0])
        writer.writeheader()
        writer.writerows(recovery_rows)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8.2))

    axes[0, 0].plot(SEARCH_PERIODS, power, color="#2f6f4f", lw=0.8)
    axes[0, 0].axhline(
        null_threshold_99,
        color="#596675",
        ls=":",
        lw=1.1,
        label="99% noise-only max-residue threshold",
    )
    axes[0, 0].axvline(
        TRUE_PERIOD_DAYS,
        color="#a8431f",
        ls="--",
        lw=1.2,
        label=f"Injected period = {TRUE_PERIOD_DAYS} d",
    )
    axes[0, 0].set(
        xlabel="Trial period [days]",
        ylabel="BLS signal residue",
        title="Period search and cadence-calibrated threshold",
    )
    axes[0, 0].legend(fontsize=7)
    axes[0, 0].grid(alpha=0.25)

    mask = np.abs(phase_days) < half_dur_days * 4
    bin_edges = np.linspace(-half_dur_days * 4, half_dur_days * 4, 60)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    binned_flux = np.array([
        flux[mask][(phase_days[mask] >= bin_edges[i]) & (phase_days[mask] < bin_edges[i + 1])].mean()
        if np.any((phase_days[mask] >= bin_edges[i]) & (phase_days[mask] < bin_edges[i + 1])) else np.nan
        for i in range(len(bin_edges) - 1)
    ])
    axes[0, 1].scatter(
        phase_days[mask] * 24,
        (flux[mask] - 1) * 1e6,
        s=3,
        color="#9fb3a8",
        alpha=0.3,
        label="Per-point flux",
    )
    axes[0, 1].plot(
        bin_centers * 24,
        (binned_flux - 1) * 1e6,
        "o-",
        color="#1f4e79",
        ms=4,
        label="Binned",
    )
    axes[0, 1].axvspan(
        -recovered_duration_hr / 2,
        recovered_duration_hr / 2,
        color="#a8431f",
        alpha=0.1,
    )
    axes[0, 1].set(
        xlabel="Hours from recovered mid-transit",
        ylabel="Relative flux [ppm]",
        title=f"Phase-folded light curve (depth SNR = {snr:.1f})",
    )
    axes[0, 1].legend(fontsize=7)
    axes[0, 1].grid(alpha=0.25)

    window_power = np.abs(
        np.exp(2j * np.pi * np.outer(1 / SEARCH_PERIODS, time)).mean(axis=1)
    ) ** 2
    axes[1, 0].plot(SEARCH_PERIODS, window_power, color="#735c8f", lw=0.8)
    axes[1, 0].set(
        xlabel="Trial period [days]",
        ylabel="Spectral-window power",
        title=f"Sampling window ({len(time) / (BASELINE_DAYS * 48):.1%} duty cycle)",
    )
    axes[1, 0].grid(alpha=0.25)

    depths = np.array([row["depth_ppm"] for row in recovery_rows])
    fractions = np.array([row["recovery_fraction"] for row in recovery_rows])
    lower = np.array([row["wilson95_low"] for row in recovery_rows])
    upper = np.array([row["wilson95_high"] for row in recovery_rows])
    axes[1, 1].errorbar(
        depths,
        fractions,
        yerr=[fractions - lower, upper - fractions],
        fmt="o-",
        color="#2f6f4f",
        capsize=3,
    )
    axes[1, 1].axvline(NOISE_PPM, color="#596675", ls=":", label="Per-point sigma")
    axes[1, 1].set(
        xscale="log",
        ylim=(-0.04, 1.04),
        xlabel="Injected transit depth [ppm]",
        ylabel="Period recovery fraction",
        title="Randomized-epoch depth sensitivity",
    )
    axes[1, 1].legend(fontsize=7)
    axes[1, 1].grid(alpha=0.25)

    fig.suptitle("Transit photometry: recovering an injected signal via Box Least Squares")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "transit_bls_recovery.png", dpi=200)

    print(f"Wrote {summary_path}")
    print(f"Wrote {sensitivity_path}")
    print(f"Wrote {FIG_DIR / 'transit_bls_recovery.png'}")
    print(f"Injected period {TRUE_PERIOD_DAYS} d -> recovered {best_period:.4f} d ({period_error_pct:.2f}% error)")
    print(f"Injected duration {TRUE_DURATION_HR} h -> recovered {recovered_duration_hr:.2f} h ({duration_error_pct:.2f}% error, searched blind over {len(DURATION_FRACTIONS)} fractions)")
    print(f"Injected depth {TRUE_DEPTH_PPM} ppm -> recovered {recovered_depth_ppm:.1f} ppm ({depth_error_pct:.2f}% error)")
    print(f"In-transit depth SNR (uncorrected for search trials): {snr:.1f}")
    print(
        f"Cadence-calibrated null: {null_exceedances}/{NULL_TRIALS} maxima "
        f"reached the observed residue; plus-one FAP={fap_plus_one:.4g}"
    )
    print(
        "Recovery fractions: "
        + ", ".join(
            f"depth={row['depth_ppm']:g}: {row['recoveries']}/{row['trials']}"
            for row in recovery_rows
        )
    )


if __name__ == "__main__":
    main()
