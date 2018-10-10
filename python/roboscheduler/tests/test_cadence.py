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
