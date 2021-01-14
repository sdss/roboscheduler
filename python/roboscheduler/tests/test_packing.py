# encoding: utf-8
#
# main.py


from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

from pytest import mark

import os
import numpy as np
import roboscheduler.cadence as cadence


class TestPacking(object):
    """Tests Packing"""

    def test_set_packing(self):
        clist = cadence.CadenceList()
        clist.reset()
        cfile = os.path.join(os.getenv('ROBOSCHEDULER_DIR'),
                             'data', 'rsc-test-cadence.fits')
        clist.fromfits(cfile)

        packing = cadence.Packing(field_cadence='faint_1x1')
        assert (packing.nepochs == 1)
        for i in np.arange(packing.nepochs, dtype=np.int32):
            assert (packing.epoch_targets[i].max() == -1)
            assert (packing.epoch_targets[i].min() == -1)
            assert (packing.epoch_nexposures[i] == 1)
            assert (packing.epoch_nused[i] == 0)
        for i in np.arange(packing.epoch_nexposures.sum(), dtype=np.int32):
            assert (packing.exposures[i] == -1)

        packing = cadence.Packing(field_cadence='faint_month_4x2')
        assert (packing.nepochs == 4)
        for i in np.arange(packing.nepochs, dtype=np.int32):
            assert (packing.epoch_targets[i].max() == -1)
            assert (packing.epoch_targets[i].min() == -1)
            assert (packing.epoch_nexposures[i] == 2)
            assert (packing.epoch_nused[i] == 0)
        for i in np.arange(packing.epoch_nexposures.sum(), dtype=np.int32):
            assert (packing.exposures[i] == -1)

    def test_check_target(self):
        clist = cadence.CadenceList()
        clist.reset()
        cfile = os.path.join(os.getenv('ROBOSCHEDULER_DIR'),
                             'data', 'rsc-test-cadence.fits')
        clist.fromfits(cfile)

        packing = cadence.Packing(field_cadence='faint_2x4')
        assert (packing.check_target(target_cadence='faint_1x1')['ok']
                is True)
        assert (packing.check_target(target_cadence='faint_month_4x2')['ok']
                is False)
        assert (packing.check_target(target_cadence='faint_2x2')['ok']
                is True)
        assert (packing.check_target(target_cadence='bright_2x2')['ok']
                is True)

        packing = cadence.Packing(field_cadence='bright_month_4x2')
        assert (packing.check_target(target_cadence='bright_1x1')['ok']
                is True)
        assert (packing.check_target(target_cadence='bright_month_4x2')['ok']
                is True)
        assert (packing.check_target(target_cadence='bright_2x2')['ok']
                is True)
        assert (packing.check_target(target_cadence='faint_2x2')['ok']
                is False)
        assert (packing.check_target(target_cadence='bright_2x4')['ok']
                is False)
        assert (packing.check_target(target_cadence='bright_month_2x2')['ok']
                is True)
        assert (packing.check_target(target_cadence='faint_month_2x2')['ok']
                is False)

    # def test_add_target(self):
    #     clist = cadence.CadenceList()
    #     clist.reset()
    #     cfile = os.path.join(os.getenv('ROBOSCHEDULER_DIR'),
    #                          'data', 'rsc-test-cadence.fits')
    #     clist.fromfits(cfile)

    #     packing = cadence.Packing(field_cadence='faint_2x4')
    #     assert(packing.add_target(target_id=1, target_cadence='faint_1x1') is
    #            True)
    #     assert(packing.add_target(target_id=1, target_cadence='faint_2x4') is
    #            False)
    #     assert(packing.epoch_targets[0][0] == 1)
    #     assert(packing.epoch_nused[0] == 1)
    #     assert(packing.epoch_nused.sum() == 1)
    #     assert(packing.exposures[0] == 1)
    #     assert(packing.exposures[1:].max() == -1)

    #     packing = cadence.Packing(field_cadence='faint_2x4')
    #     assert(packing.add_target(target_id=1, target_cadence='faint_2x4') is
    #            True)
    #     assert(packing.add_target(target_id=1, target_cadence='faint_1x4') is
    #            False)
    #     assert(packing.epoch_nused.min() == 4)
    #     assert(packing.epoch_nused.max() == 4)
    #     assert(packing.exposures.min() == 1)
    #     assert(packing.exposures.max() == 1)
