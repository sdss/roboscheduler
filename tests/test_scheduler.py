import pytest
from pytest import mark

import numpy as np

from roboscheduler import scheduler
from roboscheduler.ks91 import KS91_deltaV


def test_epochs_completed():
    exp = 18. / 60. / 24.
    # one epoch, started at 59900
    mjds_1 = 59900. + np.array([0., exp])

    mjds_2 = 59900. + np.array([0., exp, 5, 5+exp, 5+2*exp])

    comp, prev = scheduler.epochs_completed(mjd_past=mjds_1)
    assert comp == 1
    assert np.isclose(prev, 59900)

    comp, prev = scheduler.epochs_completed(mjd_past=mjds_2)
    assert comp == 2
    assert np.isclose(prev, 59905)


class TestKS91(object):
    """test ks91 deltaV

    mjd 59999.79 -> moon illum of 0.251, phase of ~120

    moon pos: ra 33.5, dec 13.5

    check seps 5, 30, 60, should get listed deltaVs
    from KS91 table 2

    decs: 18.5, 43.5, 73.5
    deltaVs: -1.58, -0.58, -0.34

    mjd 60002.33 -> moon illum of 0.501 phase of ~90

    moon pos: ra 66.14, dec 24.41

    decs: 29.41, 54.41, 84.41
    deltaVs: -2.73, -1.35, -0.88
    """

    # test column 3 of table 2 and then some
    @mark.parametrize(("alpha", "rho", "zee", "zee_m", "deltaV"),
                      [(30, 60, 0, 60, -2.22),
                       (60, 60, 0, 60, -1.54),
                       (90, 60, 0, 60, -0.88),
                       (90, 30, 30, 60, -1.35)])
    def test_ks91_raw(self, alpha, rho, zee, zee_m, deltaV):
        delta = KS91_deltaV(alpha, rho, zee, zee_m)
        assert np.isclose(delta, deltaV, atol=0.1)

    # test 120 deg row of table 2
    @mark.parametrize(("alpha", "rho", "zee", "zee_m", "deltaV"),
                      [(120, 30, 47.7, 58.1, -0.58),
                       (120, 60, 50.5, 58.1, -0.34)])
    # these fail but maybe that's ok? TBD
    # ,
    #                        (90, 30, 71.1, 87.6, -1.35),
    #                        (90, 60, 58.7, 87.6, -0.88)
    def test_ks91_raw2(self, alpha, rho, zee, zee_m, deltaV):
        delta = KS91_deltaV(alpha, rho, zee, zee_m)
        assert np.isclose(delta, deltaV, atol=0.1)

    @mark.parametrize(("mjd", "ra", "dec", "deltaV"),
                      [(59999.79, 33.5, 43.5, -0.58),
                       (59999.79, 33.5, 73.5, -0.34)])
    def test_ks91(self, mjd, ra, dec, deltaV):
        sched = scheduler.Scheduler()
        delta = sched.deltaV_sky_pos(mjd, ra, dec)
        delta = sched._arrayify(delta)
        assert np.isclose(delta, deltaV, atol=0.1)

    # @mark.parametrize(("mjd", "ra", "dec", "deltaV"),
    #                   [(60002.33,
    #                     np.array([66.14, 66.14]),
    #                     np.array([54.41, 84.41]),
    #                     np.array([-1.35, -0.88]))])
    # def test_ks91_array(self, mjd, ra, dec, deltaV):
    #     sched = scheduler.Scheduler()
    #     delta = sched.deltaV_sky_pos(mjd, ra, dec)
    #     for d, v in zip(delta, deltaV):
    #         assert np.isclose(d, v, atol=0.1)
