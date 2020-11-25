import pytest

import numpy as np
import roboscheduler.cadence as cadence


def test_add_cadence():
    clist = cadence.CadenceList()

    clist.add_cadence(name='single_1x1',
                      nepochs=1,
                      instrument='BOSS',
                      skybrightness=[1.],
                      delta=[-1.],
                      delta_min=[-1.],
                      delta_max=[-1.],
                      nexp=[1])

    clist.add_cadence(name='single_2x1',
                      nepochs=2,
                      instrument='APOGEE',
                      skybrightness=[1., 1.],
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[1, 1])

    clist.add_cadence(name='single_2x2',
                      nepochs=2,
                      instrument='APOGEE',
                      skybrightness=[1., 1.],
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[2, 2])

    assert len(clist.cadences) == 3
    assert clist.cadences['single_2x2'].nexp[1] == 2
    assert clist.cadences['single_2x1'].nexp[0] == 1
    assert clist.cadences['single_1x1'].nexp[0] == 1
    assert len(clist.cadences['single_2x1'].nexp) == 2
    assert len(clist.cadences['single_1x1'].delta_min) == 1
    assert len(clist.cadences['single_2x2'].delta_max) == 2


def test_epochs_consistency_1():
    clist = cadence.CadenceList()

    clist.add_cadence(name='timed_2x1',
                      nepochs=2,
                      skybrightness=[1., 1.],
                      instrument='BOSS',
                      delta=[0., 2.],
                      delta_min=[0., 1.],
                      delta_max=[0., 20.],
                      nexp=[1, 1])

    clist.add_cadence(name='timed_3x1',
                      nepochs=3,
                      skybrightness=[1., 1., 1.],
                      instrument='BOSS',
                      delta=[0., 2., 2.],
                      delta_min=[0., 1., 1.],
                      delta_max=[0., 20., 20.],
                      nexp=[1, 1, 1])

    tcore = clist.cadences['timed_2x1'].as_cadencecore()
    assert clist.cadences['timed_3x1'].epochs_consistency(tcore,
                                                          epochs=[0, 1]) is True
    assert clist.cadences['timed_3x1'].epochs_consistency(tcore,
                                                          epochs=[1, 2]) is True
    assert clist.cadences['timed_3x1'].epochs_consistency(tcore,
                                                          epochs=[0, 2]) is False


def test_epochs_consistency_2():
    clist = cadence.CadenceList()

    clist.add_cadence(name='single_2x1',
                      nepochs=2,
                      instrument='BOSS',
                      skybrightness=[1., 1.],
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[1, 1])

    clist.add_cadence(name='single_2x2',
                      nepochs=2,
                      instrument='BOSS',
                      skybrightness=[1., 1.],
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[2, 2])

    clist.add_cadence(name='timed_3x1',
                      nepochs=3,
                      instrument='BOSS',
                      skybrightness=[1., 1., 1.],
                      delta=[0., 2., 2.],
                      delta_min=[0., 1., 1.],
                      delta_max=[0., 20., 20.],
                      nexp=[1, 1, 1])

    tcore = clist.cadences['single_2x1'].as_cadencecore()
    assert clist.cadences['timed_3x1'].epochs_consistency(tcore,
                                                          epochs=[1, 2]) is True
    assert clist.cadences['timed_3x1'].epochs_consistency(tcore,
                                                          epochs=[0, 2]) is True
    assert clist.cadences['timed_3x1'].epochs_consistency(tcore,
                                                          epochs=[1, 2]) is True

    tcore = clist.cadences['single_2x2'].as_cadencecore()
    assert clist.cadences['timed_3x1'].epochs_consistency(tcore,
                                                          epochs=[1, 2]) is False
    assert clist.cadences['timed_3x1'].epochs_consistency(tcore,
                                                          epochs=[0, 1]) is False


def test_epochs_consistency_3():
    clist = cadence.CadenceList()

    clist.add_cadence(name='bright_2x1',
                      nepochs=2,
                      skybrightness=[1., 1.],
                      instrument='BOSS',
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[1, 1])

    clist.add_cadence(name='dark_2x1',
                      nepochs=2,
                      skybrightness=[0.35, 0.35],
                      instrument='BOSS',
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[1, 1])

    tcore = clist.cadences['bright_2x1'].as_cadencecore()
    assert clist.cadences['dark_2x1'].epochs_consistency(tcore,
                                                         epochs=[0, 1]) is True

    tcore = clist.cadences['dark_2x1'].as_cadencecore()
    assert clist.cadences['bright_2x1'].epochs_consistency(tcore,
                                                           epochs=[0, 1]) is False


def test_cadence_consistency_1():
    clist = cadence.CadenceList()

    clist.add_cadence(name='single_2x1',
                      nepochs=2,
                      skybrightness=[1., 1.],
                      instrument='BOSS',
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[1, 1])

    clist.add_cadence(name='timed_2x1',
                      nepochs=2,
                      skybrightness=[1., 1.],
                      instrument='BOSS',
                      delta=[0., 2.],
                      delta_min=[0., 1.],
                      delta_max=[0., 20.],
                      nexp=[1, 1])

    clist.add_cadence(name='timed_3x1',
                      nepochs=3,
                      skybrightness=[1., 1., 1.],
                      instrument='BOSS',
                      delta=[0., 2., 2.],
                      delta_min=[0., 1., 1.],
                      delta_max=[0., 20., 20.],
                      nexp=[1, 1, 1])

    ok, epochs_list = clist.cadence_consistency('single_2x1', 'timed_2x1')
    assert ok is True
    assert len(epochs_list) == 1
    assert epochs_list[0][0] == 0
    assert epochs_list[0][1] == 1

    ok, epochs_list = clist.cadence_consistency('timed_2x1', 'timed_3x1')
    assert ok is True
    assert len(epochs_list) == 2
    assert epochs_list[0][0] == 0
    assert epochs_list[0][1] == 1
    assert epochs_list[1][0] == 1
    assert epochs_list[1][1] == 2


def test_cadence_consistency_2():
    clist = cadence.CadenceList()

    clist.add_cadence(name='single_100x1',
                      nepochs=100,
                      instrument='BOSS',
                      skybrightness=[1.] * 100,
                      delta=[-1.] * 100,
                      delta_min=[-1.] * 100,
                      delta_max=[-1.] * 100,
                      nexp=[1] * 100)

    clist.add_cadence(name='single_1x1',
                      nepochs=1,
                      instrument='BOSS',
                      skybrightness=[1.] * 1,
                      delta=[-1.] * 1,
                      delta_min=[-1.] * 1,
                      delta_max=[-1.] * 1,
                      nexp=[1] * 1)

    clist.add_cadence(name='single_4x1',
                      nepochs=4,
                      instrument='BOSS',
                      skybrightness=[1.] * 4,
                      delta=[-1.] * 4,
                      delta_min=[-1.] * 4,
                      delta_max=[-1.] * 4,
                      nexp=[1] * 4)

    clist.add_cadence(name='single_10x1',
                      nepochs=10,
                      instrument='BOSS',
                      skybrightness=[1.] * 10,
                      delta=[-1.] * 10,
                      delta_min=[-1.] * 10,
                      delta_max=[-1.] * 10,
                      nexp=[1] * 10)

    ok, epochs_list = clist.cadence_consistency('single_1x1', 'single_100x1')
    assert ok is True
    assert len(epochs_list) == 100

    ok, epochs_list = clist.cadence_consistency('single_1x1', 'single_10x1')
    assert ok is True
    assert len(epochs_list) == 10

    ok, epochs_list = clist.cadence_consistency('single_4x1', 'single_10x1')
    assert ok is True
    assert len(epochs_list) == 210


def test_cadence_consistency_3():
    clist = cadence.CadenceList()

    clist.add_cadence(name='mwm_tess_rgb_2x1',
                      nepochs=2,
                      instrument='APOGEE',
                      skybrightness=[1., 1.],
                      delta=[0., 26.],
                      delta_min=[0., 1.],
                      delta_max=[0., 3000.],
                      nexp=[1, 1])

    clist.add_cadence(name='csc_faint_boss_1x4',
                      nepochs=1,
                      instrument='BOSS',
                      skybrightness=[0.35],
                      delta=[-1.],
                      delta_min=[-1.],
                      delta_max=[-1.],
                      nexp=[4])

    ok, epochs_list = clist.cadence_consistency('mwm_tess_rgb_2x1', 'csc_faint_boss_1x4')
    assert ok is False
