import numpy as np
import pytest
from robustness import ar1_noise, depth_gls, window_recovery_probability

def test_ar1_seed_and_stationarity():
    a=ar1_noise(100000,2,0.7,5)
    assert np.array_equal(a,ar1_noise(100000,2,0.7,5))
    assert abs(np.std(a)-2)<0.06
    assert abs(np.corrcoef(a[:-1],a[1:])[0,1]-0.7)<0.02

def test_gls_exact_box():
    mask=np.zeros(100,bool); mask[30:40]=True
    d,s=depth_gls(1-0.001*mask,mask,0.0001,0.4)
    assert abs(d-0.001)<1e-10 and s>0

def test_correlated_uncertainty():
    mask=np.zeros(300,bool);mask[50:100]=True
    _,s0=depth_gls(np.ones(300),mask,0.0001,0)
    _,s1=depth_gls(np.ones(300),mask,0.0001,0.8)
    assert s1>s0

def test_window_gaps_and_baseline():
    full=window_recovery_probability(4,0.12,60,1/48)
    gapped=window_recovery_probability(4,0.12,60,1/48,[(0,30)])
    assert full==1 and gapped<=full
    assert window_recovery_probability(4,0.12,7,1/48)==0

def test_bad_inputs():
    with pytest.raises(ValueError): ar1_noise(5,1,1)
    with pytest.raises(ValueError): window_recovery_probability(4,5,60,0.02)
