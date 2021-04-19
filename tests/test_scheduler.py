import pytest

import numpy as np

from roboscheduler import scheduler


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
