import unittest
import numpy as np
from ephemeris_diagnostics import weighted_linear_ephemeris,forecast_timing_uncertainty

class EphemerisTests(unittest.TestCase):
    def test_exact_line(self):
        e=np.arange(10.)
        p,c,r=weighted_linear_ephemeris(e,2450000+2.5*e,np.full(10,.01))
        np.testing.assert_allclose(p,[2450000,2.5],atol=1e-8)
        self.assertLess(np.max(abs(r)),1e-8)
    def test_covariance_prediction(self):
        e=np.arange(5.)
        _,c,_=weighted_linear_ephemeris(e,10+e,np.ones(5))
        self.assertAlmostEqual(float(forecast_timing_uncertainty(2,c)),1/np.sqrt(5))
        self.assertGreater(float(forecast_timing_uncertainty(100,c)),float(forecast_timing_uncertainty(2,c)))
    def test_bad_sigma(self):
        with self.assertRaises(ValueError):
            weighted_linear_ephemeris([0,1,2],[1,2,3],[1,0,1])
if __name__=="__main__": unittest.main()
