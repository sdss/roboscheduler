# encoding: utf-8
#
# main.py


from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

from pytest import mark

import os
import roboscheduler.cadence as cadence


class TestCadence(object):
    """Tests CadenceList"""

    def test_cadence_epochs(self):
        clist = cadence.CadenceList()
        clist.reset()
        clist.add_cadence(name='one', nexposures=4,
                          skybrightness=[1., 1., 1., 1.],
                          delta=[0., 0., 1., 0.],
                          delta_min=[0., 0., 0.5, 0.],
                          delta_max=[0., 0., 1.5, 0.],
                          instrument=['boss'] * 4)
        assert clist.cadences['one'].nepochs == 2
        assert clist.cadences['one'].epoch_indx[0] == 0
        assert clist.cadences['one'].epoch_indx[1] == 2

    def test_cadence_single(self):
        clist = cadence.CadenceList()
        clist.reset()
        clist.add_cadence(name='one', nexposures=1,
                          skybrightness=[1.], delta=[-1.], delta_min=[-1.],
                          delta_max=[-1.], instrument=['boss'])
        clist.add_cadence(name='two', nexposures=1,
                          skybrightness=[0.35], delta=[-1.], delta_min=[-1.],
                          delta_max=[-1.], instrument=['boss'])
        assert clist.cadence_consistency('one', 'two', return_solutions=False) == True

    def test_cadence_single_double(self):
        clist = cadence.CadenceList()
        clist.reset()
        clist.add_cadence(name='one', nexposures=1,
                          skybrightness=[1.], delta=[-1.], delta_min=[-1.],
                          delta_max=[-1.], instrument=['boss'])
        clist.add_cadence(name='two', nexposures=2,
                          skybrightness=[0.35, 0.35], delta=[0., 30.],
                          delta_min=[0., 20.],
                          delta_max=[0., 40.],
                          instrument=['boss', 'boss'])
        assert clist.cadence_consistency('one', 'two', return_solutions=False) == True
        assert clist.cadence_consistency('two', 'one', return_solutions=False) == False

    def test_cadence_single_double_skybrightness(self):
        clist = cadence.CadenceList()
        clist.reset()
        clist.add_cadence(name='bright', nexposures=1,
                          skybrightness=[1.], delta=[-1.], delta_min=[-1.],
                          delta_max=[-1.], instrument=['boss'])
        clist.add_cadence(name='dark', nexposures=1,
                          skybrightness=[0.35], delta=[-1.], delta_min=[-1.],
                          delta_max=[-1.], instrument=['boss'])
        clist.add_cadence(name='two', nexposures=2,
                          skybrightness=[1., 1.], delta=[0., 30.],
                          delta_min=[0., 20.],
                          delta_max=[0., 40.],
                          instrument=['boss', 'boss'])
        assert clist.cadence_consistency('bright', 'two',
                                         return_solutions=False) == True
        assert clist.cadence_consistency('dark', 'two',
                                         return_solutions=False) == False

    def test_cadence_double_double(self):
        clist = cadence.CadenceList()
        clist.reset()
        clist.add_cadence(name='one', nexposures=2,
                          skybrightness=[0.35, 0.35], delta=[0., 30.],
                          delta_min=[0., 19.],
                          delta_max=[0., 41.],
                          instrument=['boss', 'boss'])
        clist.add_cadence(name='two', nexposures=2,
                          skybrightness=[0.35, 0.35], delta=[0., 30.],
                          delta_min=[0., 20.],
                          delta_max=[0., 40.],
                          instrument=['boss', 'boss'])
        assert clist.cadence_consistency('one', 'two',
                                         return_solutions=False) is True
        assert clist.cadence_consistency('two', 'one',
                                         return_solutions=False) is False

    def test_cadence_file(self):
        clist = cadence.CadenceList()
        clist.reset()
        cfile = os.path.join(os.getenv('ROBOSCHEDULER_DIR'),
                             'data', 'rsc-test-cadence.fits')
        clist.fromfits(cfile)

        assert clist.ncadences == 12

        assert (clist.cadence_consistency('bright_1x1', 'faint_1x1',
                                          return_solutions=False) is True)

        assert (clist.cadence_consistency('faint_1x1', 'bright_1x1',
                                          return_solutions=False) is False)

        assert (clist.cadence_consistency('faint_1x1', 'faint_1x4',
                                          return_solutions=False) is True)

        assert (clist.cadence_consistency('faint_1x4', 'faint_1x1',
                                          return_solutions=False) is False)

        assert (clist.cadence_consistency('faint_1x4', 'faint_1x4',
                                          return_solutions=False) is True)

        assert (clist.cadence_consistency('faint_1x4', 'faint_1x4',
                                          return_solutions=False) is True)

        assert (clist.cadence_consistency('faint_2x2', 'faint_month_2x2',
                                          return_solutions=False) is True)

        assert (clist.cadence_consistency('faint_2x2', 'faint_month_4x2',
                                          return_solutions=False) is True)

        assert (clist.cadence_consistency('faint_month_2x2', 'faint_month_4x2',
                                          return_solutions=False) is True)

        assert (clist.cadence_consistency('faint_month_2x2', 'faint_month_4x2',
                                          return_solutions=False) is True)

        assert (clist.cadence_consistency('faint_1x4', 'faint_month_4x2',
                                          return_solutions=False) is False)
