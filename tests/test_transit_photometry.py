"""Executable checks on the transit-model and BLS-search physics, not
just the illustrative injection-recovery run in the demo script itself.
"""

import numpy as np
import transit_photometry_demo as tpd


def test_transit_model_depth_at_midtransit():
    period, t0, depth_ppm, duration_hr = 4.0, 2.0, 500.0, 3.0
    flux = tpd.transit_model(np.array([t0]), period, t0, depth_ppm, duration_hr)
    assert np.isclose(flux[0], 1.0 - depth_ppm * 1e-6, atol=1e-9)


def test_cadence_is_deterministic_gapped_and_ordered():
    first = tpd.make_cadence(np.random.default_rng(11))
    second = tpd.make_cadence(np.random.default_rng(11))
    full_count = int(tpd.BASELINE_DAYS * 24 * 60 / tpd.CADENCE_MIN)
    assert np.array_equal(first, second)
    assert np.all(np.diff(first) > 0)
    assert 0.9 < len(first) / full_count < 0.96
    assert not np.any((first >= 14.5) & (first < 15.25))


def test_wilson_interval_handles_boundary_counts():
    low_zero, high_zero = tpd.wilson_interval(0, 20)
    low_full, high_full = tpd.wilson_interval(20, 20)
    assert abs(low_zero) < 1e-12
    assert 0.15 < high_zero < 0.17
    assert 0.83 < low_full < 0.85
    assert abs(high_full - 1) < 1e-12


def test_transit_model_flat_out_of_transit():
    period, t0, depth_ppm, duration_hr = 4.0, 2.0, 500.0, 3.0
    # A full period away from t0, well outside any transit or taper window.
    far_time = np.array([t0 + period / 2])
    flux = tpd.transit_model(far_time, period, t0, depth_ppm, duration_hr)
    assert np.isclose(flux[0], 1.0, atol=1e-9)


def test_transit_model_conserves_baseline_flux_scale():
    # Every point is either exactly 1.0 (out of transit) or below 1.0
    # (in transit / ingress-egress) -- never above baseline.
    time = np.linspace(0, 20, 5000)
    flux = tpd.transit_model(time, period=4.35, t0=1.7, depth_ppm=850.0, duration_hr=2.8)
    assert flux.max() <= 1.0 + 1e-12
    assert flux.min() >= 1.0 - 850e-6 - 1e-9


def test_best_box_recovers_injected_depth_noiseless():
    period, t0, depth_ppm, duration_hr = 4.35, 1.7, 850.0, 2.8
    time = np.arange(0, 60.0, 30.0 / (24 * 60))
    flux = tpd.transit_model(time, period, t0, depth_ppm, duration_hr)
    result = tpd.best_box_for_period(time, flux, period, duration_hr)
    # Noiseless injection at the true period: recovered depth should match
    # the injected depth almost exactly (ingress/egress smoothing is the
    # only source of small disagreement).
    assert abs(result["depth"] * 1e6 - depth_ppm) / depth_ppm < 0.05


def test_box_least_squares_recovers_period_across_seeds():
    period, t0, depth_ppm, duration_hr = 4.35, 1.7, 850.0, 2.8
    time = np.arange(0, 60.0, 30.0 / (24 * 60))
    true_flux = tpd.transit_model(time, period, t0, depth_ppm, duration_hr)
    period_grid = np.arange(3.5, 5.5, 0.01)

    for seed in (1, 2, 3):
        rng = np.random.default_rng(seed)
        flux = true_flux + rng.normal(0, 100e-6, size=time.size)
        power = tpd.box_least_squares(time, flux, period_grid, tpd.DURATION_FRACTIONS)
        recovered_period = period_grid[np.argmax(power)]
        assert abs(recovered_period - period) < 0.05, f"seed {seed}: recovered {recovered_period}, injected {period}"


def test_kepler_class_noise_regime_is_realistic():
    # Sanity check that the noise level used in the demo matches Kepler's
    # own published per-30-min photometric precision order of magnitude
    # for a quiet Sun-like star, not an arbitrarily easy/hard value.
    assert 50.0 <= tpd.NOISE_PPM <= 200.0
