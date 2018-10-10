# encoding: utf-8
#
# main.py


from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

from pytest import mark

import roboscheduler.cadence as cadence


class TestCadence(object):
    """Tests CadenceList"""

    def test_cadence_epochs(self):
        clist = cadence.CadenceList()
        clist.reset()
        clist.add_cadence(name='one', nexposures=4,
                          lunation=[1., 1., 1., 1.],
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
                          lunation=[1.], delta=[-1.], delta_min=[-1.],
                          delta_max=[-1.], instrument=['boss'])
        clist.add_cadence(name='two', nexposures=1,
                          lunation=[0.35], delta=[-1.], delta_min=[-1.],
                          delta_max=[-1.], instrument=['boss'])
        assert clist.cadence_consistency('one', 'two', return_solutions=False) == True

    def test_cadence_single_double(self):
        clist = cadence.CadenceList()
        clist.reset()
        clist.add_cadence(name='one', nexposures=1,
                          lunation=[1.], delta=[-1.], delta_min=[-1.],
                          delta_max=[-1.], instrument=['boss'])
        clist.add_cadence(name='two', nexposures=2,
                          lunation=[0.35, 0.35], delta=[0., 30.],
                          delta_min=[0., 20.],
                          delta_max=[0., 40.],
                          instrument=['boss', 'boss'])
        assert clist.cadence_consistency('one', 'two', return_solutions=False) == True
        assert clist.cadence_consistency('two', 'one', return_solutions=False) == False

    def test_cadence_single_double_lunation(self):
        clist = cadence.CadenceList()
        clist.reset()
        clist.add_cadence(name='bright', nexposures=1,
                          lunation=[1.], delta=[-1.], delta_min=[-1.],
                          delta_max=[-1.], instrument=['boss'])
        clist.add_cadence(name='dark', nexposures=1,
                          lunation=[0.35], delta=[-1.], delta_min=[-1.],
                          delta_max=[-1.], instrument=['boss'])
        clist.add_cadence(name='two', nexposures=2,
                          lunation=[1., 1.], delta=[0., 30.],
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
                          lunation=[0.35, 0.35], delta=[0., 30.],
                          delta_min=[0., 19.],
                          delta_max=[0., 41.],
                          instrument=['boss', 'boss'])
        clist.add_cadence(name='two', nexposures=2,
                          lunation=[0.35, 0.35], delta=[0., 30.],
                          delta_min=[0., 20.],
                          delta_max=[0., 40.],
                          instrument=['boss', 'boss'])
        assert clist.cadence_consistency('one', 'two', return_solutions=False) == True
        assert clist.cadence_consistency('two', 'one', return_solutions=False) == False

    def test_cadence_pack(self):
        clist = cadence.CadenceList()
        clist.reset()
        clist.add_cadence(name='one', nexposures=1,
                          lunation=[1.], delta=[-1.], delta_min=[-1.],
                          delta_max=[-1.], instrument=['boss'])
        clist.add_cadence(name='two', nexposures=2,
                          lunation=[1., 1.],
                          delta=[0., 1.],
                          delta_min=[0., 0.4],
                          delta_max=[0., 1.6], instrument=['boss'])
        clist.add_cadence(name='three', nexposures=4,
                          lunation=[1., 1., 1., 1.],
                          delta=[0., 0., 1., 0.],
                          delta_min=[0., 0., 0.5, 0.],
                          delta_max=[0., 0., 1.5, 0.],
                          instrument=['boss'] * 4)

        (epoch_targets, exposure_targets) = (
            clist.pack_targets(target_cadences=['one', 'two', 'one', 'two'],
                               value=[2, 2, 1, 1],
                               field_cadence='three'))

        assert epoch_targets[0][0] == 1
        assert epoch_targets[0][1] == 2
        assert epoch_targets[1][0] == 0
        assert epoch_targets[1][1] == 1
        assert list(exposure_targets) == [1, 2, 0, 1]

        (epoch_targets, exposure_targets) = (
            clist.pack_targets(target_cadences=['one', 'two', 'one', 'two'],
                               value=[2, 2, 1, 10],
                               field_cadence='three'))

        assert epoch_targets[0][0] == 2
        assert epoch_targets[0][1] == 3
        assert epoch_targets[1][0] == 0
        assert epoch_targets[1][1] == 3
        assert list(exposure_targets) == [2, 3, 0, 3]

    def test_cadence_pack_greedy(self):
        clist = cadence.CadenceList()
        clist.reset()
        clist.add_cadence(name='one', nexposures=1,
                          lunation=[1.], delta=[-1.], delta_min=[-1.],
                          delta_max=[-1.], instrument=['boss'])
        clist.add_cadence(name='two', nexposures=2,
                          lunation=[1., 1.],
                          delta=[0., 1.],
                          delta_min=[0., 0.4],
                          delta_max=[0., 1.6], instrument=['boss'])
        clist.add_cadence(name='three', nexposures=4,
                          lunation=[1., 1., 1., 1.],
                          delta=[0., 0., 1., 0.],
                          delta_min=[0., 0., 0.5, 0.],
                          delta_max=[0., 0., 1.5, 0.],
                          instrument=['boss'] * 4)

        (epoch_targets, exposure_targets) = (
            clist.pack_targets_greedy(target_cadences=['one', 'two',
                                                       'one', 'two'],
                                      value=[2, 2, 1, 1],
                                      field_cadence='three'))

        assert epoch_targets[0][0] == 1
        assert epoch_targets[0][1] == 0
        assert epoch_targets[1][0] == 1
        assert epoch_targets[1][1] == 2
        assert list(exposure_targets) == [1, 0, 1, 2]

        (epoch_targets, exposure_targets) = (
            clist.pack_targets_greedy(target_cadences=['one', 'two',
                                                       'one', 'two'],
                                      value=[2, 2, 1, 10],
                                      field_cadence='three'))

        assert epoch_targets[0][0] == 3
        assert epoch_targets[0][1] == 1
        assert epoch_targets[1][0] == 3
        assert epoch_targets[1][1] == 1
        assert list(exposure_targets) == [3, 1, 3, 1]
