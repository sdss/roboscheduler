import os
import pytest

import numpy as np
import roboscheduler.cadence as cadence


def add_cadence_single_nxm(n=1, m=1):
    clist = cadence.CadenceList()
    clist.add_cadence(name='single_{n}x{m}'.format(n=n, m=m),
                      nepochs=n,
                      skybrightness=[1.] * n,
                      delta=[-1.] * n,
                      delta_min=[-1.] * n,
                      delta_max=[-1.] * n,
                      nexp=[m] * n,
                      max_length=[0.] * n,
                      min_moon_sep=[15.] * n,
                      min_deltav_ks91=[-2.5] * n,
                      min_twilight_ang=[8] * n,
                      max_airmass=[2.] * n,
                      obsmode_pk=['bright_time'] * n)
    return


def add_cadence_mixedreverse2_nxm(n=2, m=1):
    clist = cadence.CadenceList()
    clist.add_cadence(name='mixedreverse2_{n}x{m}'.format(n=n, m=m),
                      nepochs=n,
                      skybrightness=[1.00, 1.00] + [0.35] * (n - 2),
                      delta=[-1.] * n,
                      delta_min=[-1.] * n,
                      delta_max=[-1.] * n,
                      nexp=[m] * n,
                      max_length=[0.] * n,
                      min_moon_sep=[15.] * n,
                      min_deltav_ks91=[-2.5] * 2 + [-1.5] * (n - 2),
                      min_twilight_ang=[8] * 2 + [15.] * (n - 2),
                      max_airmass=[2.] * 2 + [1.4] * (n - 2),
                      obsmode_pk=['bright_time'] * 2 + ['dark_plane'] * (n - 2))
    return


def add_cadence_mixed2_nxm(n=2, m=1):
    clist = cadence.CadenceList()
    clist.add_cadence(name='mixed2_{n}x{m}'.format(n=n, m=m),
                      nepochs=n,
                      skybrightness=[0.35, 0.35] + [1.] * (n - 2),
                      delta=[0., 3.] + [-1.] +  [30.] * (n - 3),
                      delta_min=[0., 0.5] + [-1.] +  [0.5] * (n - 3),
                      delta_max=[0., 1800.] + [-1.] +  [1800.] * (n - 3),
                      nexp=[m] * n,
                      max_length=[0.] * n,
                      min_moon_sep=[15.] * n,
                      min_deltav_ks91=[-1.5] * 2 + [-2.5] * (n - 2),
                      min_twilight_ang=[15] * 2 + [8.] * (n - 2),
                      max_airmass=[1.4] * 2 + [2.] * (n - 2),
                      obsmode_pk=['dark_plane'] * 2 + ['bright_time'] * (n - 2))
    return


def test_add_cadence():
    clist = cadence.CadenceList()
    clist.reset()

    clist.add_cadence(name='single_1x1',
                      nepochs=1,
                      skybrightness=[1.],
                      delta=[-1.],
                      delta_min=[-1.],
                      delta_max=[-1.],
                      nexp=[1],
                      max_length=[1.],
                      min_deltav_ks91=[-2.5],
                      max_airmass=[2],
                      min_twilight_ang=[8],
                      min_moon_sep=[15])

    clist.add_cadence(name='single_2x1',
                      nepochs=2,
                      skybrightness=[1., 1.],
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[1, 1],
                      max_length=[1., 1.],
                      min_deltav_ks91=[-2.5, -2.5],
                      max_airmass=[2, 2],
                      min_twilight_ang=[8, 8],
                      min_moon_sep=[15, 15])

    clist.add_cadence(name='single_2x2',
                      nepochs=2,
                      skybrightness=[1., 1.],
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[2, 2],
                      max_length=[1., 1.],
                      min_deltav_ks91=[-2.5, -2.5],
                      max_airmass=[2, 2],
                      min_twilight_ang=[8, 8],
                      min_moon_sep=[15, 15])

    assert len(clist.cadences) == 3
    assert clist.cadences['single_2x2'].nexp[1] == 2
    assert clist.cadences['single_2x1'].nexp[0] == 1
    assert clist.cadences['single_1x1'].nexp[0] == 1
    assert len(clist.cadences['single_2x1'].nexp) == 2
    assert len(clist.cadences['single_1x1'].delta_min) == 1
    assert len(clist.cadences['single_2x2'].delta_max) == 2


def test_add_cadence_cfg():
    clist = cadence.CadenceList()
    clist.reset()

    test_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        'test_cfginput.cfg',
        )

    clist.fromcfg(test_file)

    assert clist.cadences['single_1x1'].nepochs == 1
    assert clist.cadences['single_1x1'].nexp[0] == 1
    assert clist.cadences['single_1x1'].delta[0] == -1.
    assert clist.cadences['single_1x1'].delta_min[0] == -1.
    assert clist.cadences['single_1x1'].delta_max[0] == -1.
    assert np.abs(clist.cadences['single_1x1'].skybrightness[0] - 0.35) < 0.001
    assert clist.cadences['single_1x1'].max_length[0] == 100.

    assert clist.cadences['single_3xY'].nepochs == 3
    assert clist.cadences['single_3xY'].nexp[0] == 2
    assert clist.cadences['single_3xY'].nexp[1] == 3
    assert clist.cadences['single_3xY'].nexp[2] == 4
    assert clist.cadences['single_3xY'].delta[0] == 3.
    assert clist.cadences['single_3xY'].delta[1] == 4.
    assert clist.cadences['single_3xY'].delta[2] == 5.
    assert clist.cadences['single_3xY'].delta_min[0] == 6.
    assert clist.cadences['single_3xY'].delta_min[1] == 7.
    assert clist.cadences['single_3xY'].delta_min[2] == 8.
    assert clist.cadences['single_3xY'].delta_max[0] == 9.
    assert clist.cadences['single_3xY'].delta_max[1] == 10.
    assert clist.cadences['single_3xY'].delta_max[2] == 11.
    assert clist.cadences['single_3xY'].skybrightness[0] == 1.
    assert clist.cadences['single_3xY'].skybrightness[1] == 2.
    assert clist.cadences['single_3xY'].skybrightness[2] == 3.
    assert clist.cadences['single_3xY'].max_length[0] == 10.
    assert clist.cadences['single_3xY'].max_length[1] == 100.
    assert clist.cadences['single_3xY'].max_length[2] == 1000.

    arr = clist.toarray()
    assert len(arr) == 2
    assert arr['NEPOCHS'][0] == 1
    return


def test_epochs_consistency_1():
    clist = cadence.CadenceList()
    clist.reset()

    clist.add_cadence(name='timed_2x1',
                      nepochs=2,
                      skybrightness=[1., 1.],
                      delta=[0., 2.],
                      delta_min=[0., 1.],
                      delta_max=[0., 20.],
                      nexp=[1, 1],
                      max_length=[1., 1.],
                      min_deltav_ks91=[-2.5, -2.5],
                      max_airmass=[2, 2],
                      min_twilight_ang=[8, 8],
                      min_moon_sep=[15, 15])

    clist.add_cadence(name='timed_3x1',
                      nepochs=3,
                      skybrightness=[1., 1., 1.],
                      delta=[0., 2., 2.],
                      delta_min=[0., 1., 1.],
                      delta_max=[0., 20., 20.],
                      nexp=[1, 1, 1],
                      max_length=[1., 1., 1.],
                      min_deltav_ks91=[-2.5, -2.5, -2.5],
                      max_airmass=[2, 2, 2],
                      min_twilight_ang=[8, 8, 8],
                      min_moon_sep=[15, 15, 15])

    tcore = clist.cadences['timed_2x1'].as_cadencecore()
    assert clist.cadences['timed_3x1'].epochs_consistency(tcore,
                                                          epochs=[0, 1],
                                                          skybrightness_only=False, inorder=0) is True
    assert clist.cadences['timed_3x1'].epochs_consistency(tcore,
                                                          epochs=[1, 2],
                                                          skybrightness_only=False, inorder=0) is True
    assert clist.cadences['timed_3x1'].epochs_consistency(tcore,
                                                          epochs=[0, 2],
                                                          skybrightness_only=False, inorder=0) is False


def test_epochs_consistency_2():
    clist = cadence.CadenceList()
    clist.reset()

    clist.add_cadence(name='single_2x1',
                      nepochs=2,
                      skybrightness=[1., 1.],
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[1, 1],
                      max_length=[1., 1.],
                      min_deltav_ks91=[-2.5, -2.5],
                      max_airmass=[2, 2],
                      min_twilight_ang=[8, 8],
                      min_moon_sep=[15, 15])

    clist.add_cadence(name='single_2x2',
                      nepochs=2,
                      skybrightness=[1., 1.],
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[2, 2],
                      max_length=[1., 1.],
                      min_deltav_ks91=[-2.5, -2.5],
                      max_airmass=[2, 2],
                      min_twilight_ang=[8, 8],
                      min_moon_sep=[15, 15])

    clist.add_cadence(name='timed_3x1',
                      nepochs=3,
                      skybrightness=[1., 1., 1.],
                      delta=[0., 2., 2.],
                      delta_min=[0., 1., 1.],
                      delta_max=[0., 20., 20.],
                      nexp=[1, 1, 1],
                      max_length=[1., 1., 1.],
                      min_deltav_ks91=[-2.5, -2.5, -2.5],
                      max_airmass=[2, 2, 2],
                      min_twilight_ang=[8, 8, 8],
                      min_moon_sep=[15, 15, 15])

    tcore = clist.cadences['single_2x1'].as_cadencecore()
    assert clist.cadences['timed_3x1'].epochs_consistency(tcore,
                                                          epochs=[1, 2],
                                                          skybrightness_only=False, inorder=0) is True
    assert clist.cadences['timed_3x1'].epochs_consistency(tcore,
                                                          epochs=[0, 2],
                                                          skybrightness_only=False, inorder=0) is True
    assert clist.cadences['timed_3x1'].epochs_consistency(tcore,
                                                          epochs=[1, 2],
                                                          skybrightness_only=False, inorder=0) is True

    tcore = clist.cadences['single_2x2'].as_cadencecore()
    assert clist.cadences['timed_3x1'].epochs_consistency(tcore,
                                                          epochs=[1, 2],
                                                          skybrightness_only=False, inorder=0) is False
    assert clist.cadences['timed_3x1'].epochs_consistency(tcore,
                                                          epochs=[0, 1],
                                                          skybrightness_only=False, inorder=0) is False
    return


def test_epochs_consistency_3():
    clist = cadence.CadenceList()
    clist.reset()

    clist.add_cadence(name='bright_2x1',
                      nepochs=2,
                      skybrightness=[1., 1.],
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[1, 1],
                      max_length=[1., 1.],
                      min_deltav_ks91=[-2.5, -2.5],
                      max_airmass=[2, 2],
                      min_twilight_ang=[8, 8],
                      min_moon_sep=[15, 15])

    clist.add_cadence(name='dark_2x1',
                      nepochs=2,
                      skybrightness=[0.35, 0.35],
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[1, 1],
                      max_length=[1., 1.],
                      min_deltav_ks91=[-1.5, -1.5],
                      max_airmass=[1.4, 1.4],
                      min_twilight_ang=[15, 15],
                      min_moon_sep=[35, 35])

    tcore = clist.cadences['bright_2x1'].as_cadencecore()
    assert clist.cadences['dark_2x1'].epochs_consistency(tcore,
                                                         epochs=[0, 1],
                                                         skybrightness_only=False, inorder=0) is True

    tcore = clist.cadences['dark_2x1'].as_cadencecore()
    assert clist.cadences['bright_2x1'].epochs_consistency(tcore,
                                                           epochs=[0, 1],
                                                           skybrightness_only=False, inorder=0) is False


def test_epochs_consistency_4():
    clist = cadence.CadenceList()
    clist.reset()

    add_cadence_single_nxm(n=4, m=1)
    add_cadence_single_nxm(n=2, m=1)
    add_cadence_mixed2_nxm(n=5, m=1)
    add_cadence_mixed2_nxm(n=8, m=1)
    add_cadence_mixed2_nxm(n=15, m=1)
    add_cadence_mixedreverse2_nxm(n=5, m=1)

    tcore = clist.cadences['single_4x1'].as_cadencecore()
    assert clist.cadences['single_4x1'].epochs_consistency(tcore,
                                                           epochs=[0, 1, 2, 3],
                                                           skybrightness_only=False, inorder=1) is True


    assert clist.cadences['single_4x1'].epochs_consistency(tcore,
                                                           epochs=[0, 1, 2],
                                                           skybrightness_only=False, inorder=1) is True


    tcore = clist.cadences['single_2x1'].as_cadencecore()
    assert clist.cadences['single_4x1'].epochs_consistency(tcore,
                                                           epochs=[0, 1],
                                                           skybrightness_only=False, inorder=1) is True

    tcore = clist.cadences['single_2x1'].as_cadencecore()
    assert clist.cadences['single_4x1'].epochs_consistency(tcore,
                                                           epochs=[2, 3],
                                                           skybrightness_only=False, inorder=0) is True

    tcore = clist.cadences['mixed2_5x1'].as_cadencecore()
    assert clist.cadences['mixed2_5x1'].epochs_consistency(tcore,
                                                           epochs=[1],
                                                           skybrightness_only=False, inorder=1) is True

    tcore = clist.cadences['mixed2_5x1'].as_cadencecore()
    assert clist.cadences['mixed2_5x1'].epochs_consistency(tcore,
                                                           epochs=[0, 2, 3],
                                                           skybrightness_only=False, inorder=1) is True

    tcore = clist.cadences['mixed2_8x1'].as_cadencecore()
    assert clist.cadences['mixed2_5x1'].epochs_consistency(tcore,
                                                           epochs=[0, 2, 3],
                                                           skybrightness_only=False, inorder=1) is True

    tcore = clist.cadences['mixed2_8x1'].as_cadencecore()
    assert clist.cadences['mixed2_15x1'].epochs_consistency(tcore,
                                                            epochs=[0, 2, 3],
                                                            skybrightness_only=False, inorder=1) is True

    tcore = clist.cadences['mixedreverse2_5x1'].as_cadencecore()
    assert clist.cadences['mixedreverse2_5x1'].epochs_consistency(tcore,
                                                                  epochs=[0, 2],
                                                                  skybrightness_only=False, inorder=1) is True

    return


def test_skybrightness_only():
    clist = cadence.CadenceList()
    clist.reset()

    clist.add_cadence(name='dark_2x1',
                      nepochs=2,
                      skybrightness=[0.35, 0.35],
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[1, 1],
                      max_length=[1., 1.],
                      min_deltav_ks91=[-1.5, -1.5],
                      max_airmass=[2, 2],
                      min_twilight_ang=[8, 8],
                      min_moon_sep=[15, 15])

    clist.add_cadence(name='darker_2x1',
                      nepochs=2,
                      skybrightness=[0.35, 0.35],
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[1, 1],
                      max_length=[1., 1.],
                      min_deltav_ks91=[-1.5, -1.5],
                      max_airmass=[1.4, 1.4],
                      min_twilight_ang=[15, 15],
                      min_moon_sep=[35, 35])

    tcore = clist.cadences['darker_2x1'].as_cadencecore()
    assert clist.cadences['dark_2x1'].epochs_consistency(tcore,
                                                         epochs=[0, 1],
                                                         skybrightness_only=False, inorder=0) is False

    tcore = clist.cadences['darker_2x1'].as_cadencecore()
    assert clist.cadences['dark_2x1'].epochs_consistency(tcore,
                                                         epochs=[0, 1],
                                                         skybrightness_only=True, inorder=0) is True

    return

def test_mixed_consistency():
    clist = cadence.CadenceList()
    clist.reset()

    clist.add_cadence(name='mixed_2x1',
                      nepochs=2,
                      skybrightness=[0.35, 1.],
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[1, 1],
                      max_length=[1., 1.],
                      min_deltav_ks91=[-1.5, -2.5],
                      max_airmass=[1.4, 2],
                      min_twilight_ang=[15, 8],
                      min_moon_sep=[15, 15])

    clist.add_cadence(name='mixedup_2x1',
                      nepochs=2,
                      skybrightness=[1.0, 0.35],
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[1, 1],
                      max_length=[1., 1.],
                      min_deltav_ks91=[-2.5, -1.5],
                      max_airmass=[2., 1.4],
                      min_twilight_ang=[8, 15],
                      min_moon_sep=[15, 15])

    clist.add_cadence(name='dark_2x1',
                      nepochs=2,
                      skybrightness=[0.35, 0.35],
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[1, 1],
                      max_length=[1., 1.],
                      min_deltav_ks91=[-1.5, -1.5],
                      max_airmass=[1.4, 1.4],
                      min_twilight_ang=[15, 15],
                      min_moon_sep=[35, 35])

    clist.add_cadence(name='bright_2x1',
                      nepochs=2,
                      skybrightness=[1.00, 1.00],
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[1, 1],
                      max_length=[1., 1.],
                      min_deltav_ks91=[-2.5, -2.5],
                      max_airmass=[2., 2.],
                      min_twilight_ang=[8, 8],
                      min_moon_sep=[35, 35])

    tcore = clist.cadences['bright_2x1'].as_cadencecore()
    assert clist.cadences['dark_2x1'].epochs_consistency(tcore,
                                                         epochs=[0, 1],
                                                         skybrightness_only=False, inorder=0) is True

    tcore = clist.cadences['dark_2x1'].as_cadencecore()
    assert clist.cadences['bright_2x1'].epochs_consistency(tcore,
                                                           epochs=[0, 1],
                                                           skybrightness_only=False, inorder=0) is False

    tcore = clist.cadences['mixed_2x1'].as_cadencecore()
    assert clist.cadences['bright_2x1'].epochs_consistency(tcore,
                                                           epochs=[0, 1],
                                                           skybrightness_only=False, inorder=0) is False

    tcore = clist.cadences['bright_2x1'].as_cadencecore()
    assert clist.cadences['mixed_2x1'].epochs_consistency(tcore,
                                                          epochs=[0, 1],
                                                          skybrightness_only=False, inorder=0) is True

    tcore = clist.cadences['dark_2x1'].as_cadencecore()
    assert clist.cadences['mixed_2x1'].epochs_consistency(tcore,
                                                          epochs=[0, 1],
                                                          skybrightness_only=False, inorder=0) is False

    tcore = clist.cadences['mixedup_2x1'].as_cadencecore()
    assert clist.cadences['mixed_2x1'].epochs_consistency(tcore,
                                                          epochs=[0, 1],
                                                          skybrightness_only=False, inorder=0) is False

    return


def test_epochs_consistency_5():
    clist = cadence.CadenceList()
    clist.reset()

    clist.add_cadence(name='bright_3x2',
                      nepochs=3,
                      skybrightness=[1., 1., 1.],
                      delta=[-1., -1., -1.],
                      delta_min=[-1., -1., -1.],
                      delta_max=[-1., -1., -1.],
                      nexp=[2, 2, 2],
                      max_length=[1., 1., 1.],
                      min_deltav_ks91=[-2.5, -2.5, -2.5],
                      max_airmass=[2, 2, 2],
                      min_twilight_ang=[8, 8, 8],
                      min_moon_sep=[15, 15, 15])

    clist.add_cadence(name='bright_3x1',
                      nepochs=3,
                      skybrightness=[1., 1., 1.],
                      delta=[-1., -1., -1.],
                      delta_min=[-1., -1., -1.],
                      delta_max=[-1., -1., -1.],
                      nexp=[1, 1, 1],
                      max_length=[1., 1., 1.],
                      min_deltav_ks91=[-2.5, -2.5, -2.5],
                      max_airmass=[2, 2, 2],
                      min_twilight_ang=[8, 8, 8],
                      min_moon_sep=[15, 15, 15])

    tcore = clist.cadences['bright_3x1'].as_cadencecore()
    assert clist.cadences['bright_3x2'].epochs_consistency(tcore,
                                                           epochs=[0, 1, 2],
                                                           skybrightness_only=False, inorder=0) is True

    tcore = clist.cadences['bright_3x1'].as_cadencecore()
    assert clist.cadences['bright_3x2'].epochs_consistency(tcore,
                                                           epochs=[0, 0, 2],
                                                           skybrightness_only=False, inorder=0) is True

    tcore = clist.cadences['bright_3x1'].as_cadencecore()
    assert clist.cadences['bright_3x2'].epochs_consistency(tcore,
                                                           epochs=[0, 1, 1],
                                                           skybrightness_only=False, inorder=0) is True

    tcore = clist.cadences['bright_3x1'].as_cadencecore()
    assert clist.cadences['bright_3x2'].epochs_consistency(tcore,
                                                           epochs=[0, 0, 0],
                                                           skybrightness_only=False, inorder=0) is False


def test_exposure_consistency():
    clist = cadence.CadenceList()
    clist.reset()

    clist.add_cadence(name='bright_3x2',
                      nepochs=3,
                      skybrightness=[1., 1., 1.],
                      delta=[-1., -1., -1.],
                      delta_min=[-1., -1., -1.],
                      delta_max=[-1., -1., -1.],
                      nexp=[2, 2, 2],
                      max_length=[1., 1., 1.],
                      min_deltav_ks91=[-2.5, -2.5, -2.5],
                      max_airmass=[2, 2, 2],
                      min_twilight_ang=[8, 8, 8],
                      min_moon_sep=[15, 15, 15])

    clist.add_cadence(name='bright_3x1',
                      nepochs=3,
                      skybrightness=[1., 1., 1.],
                      delta=[-1., -1., -1.],
                      delta_min=[-1., -1., -1.],
                      delta_max=[-1., -1., -1.],
                      nexp=[1, 1, 1],
                      max_length=[1., 1., 1.],
                      min_deltav_ks91=[-2.5, -2.5, -2.5],
                      max_airmass=[2, 2, 2],
                      min_twilight_ang=[8, 8, 8],
                      min_moon_sep=[15, 15, 15])

    assert clist.exposure_consistency('bright_3x1', 'bright_3x2',
                                      [0, 1, 2]) is True

    assert clist.exposure_consistency('bright_3x1', 'bright_3x2',
                                      [0, 0, 2]) is True

    assert clist.exposure_consistency('bright_3x1', 'bright_3x2',
                                      [0, 0, 0]) is False

    assert clist.exposure_consistency('bright_3x1', 'bright_3x2',
                                      [0, 2, 2]) is True

    assert clist.exposure_consistency('bright_3x1', 'bright_3x2',
                                      [2, 2, 2]) is False

    assert clist.exposure_consistency('bright_3x2', 'bright_3x2',
                                      [0, 1, 2, 3, 4]) is False

    assert clist.exposure_consistency('bright_3x2', 'bright_3x2',
                                      [0, 1, 2, 3, 4, 5]) is True


def test_cadence_consistency_1():
    clist = cadence.CadenceList()
    clist.reset()

    clist.add_cadence(name='single_2x1',
                      nepochs=2,
                      skybrightness=[1., 1.],
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[1, 1],
                      max_length=[1., 1.],
                      min_deltav_ks91=[-2.5, -2.5],
                      max_airmass=[2, 2],
                      min_twilight_ang=[8, 8],
                      min_moon_sep=[15, 15])

    clist.add_cadence(name='timed_2x1',
                      nepochs=2,
                      skybrightness=[1., 1.],
                      delta=[0., 2.],
                      delta_min=[0., 1.],
                      delta_max=[0., 20.],
                      nexp=[1, 1],
                      max_length=[1., 1.],
                      min_deltav_ks91=[-2.5, -2.5],
                      max_airmass=[2, 2],
                      min_twilight_ang=[8, 8],
                      min_moon_sep=[15, 15])

    clist.add_cadence(name='timed_3x1',
                      nepochs=3,
                      skybrightness=[1., 1., 1.],
                      delta=[0., 2., 2.],
                      delta_min=[0., 1., 1.],
                      delta_max=[0., 20., 20.],
                      nexp=[1, 1, 1],
                      max_length=[1., 1., 1.],
                      min_deltav_ks91=[-2.5, -2.5, -2.5],
                      max_airmass=[2, 2, 2],
                      min_twilight_ang=[8, 8, 8],
                      min_moon_sep=[15, 15, 15])

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
    clist.reset()

    clist.add_cadence(name='single_100x1',
                      nepochs=100,
                      skybrightness=[1.] * 100,
                      delta=[-1.] * 100,
                      delta_min=[-1.] * 100,
                      delta_max=[-1.] * 100,
                      nexp=[1] * 100,
                      max_length=[1.] * 100,
                      min_deltav_ks91=[-2.5] * 100,
                      max_airmass=[2] * 100,
                      min_twilight_ang=[8] * 100,
                      min_moon_sep=[15] * 100)

    clist.add_cadence(name='single_1x1',
                      nepochs=1,
                      skybrightness=[1.] * 1,
                      delta=[-1.] * 1,
                      delta_min=[-1.] * 1,
                      delta_max=[-1.] * 1,
                      nexp=[1] * 1,
                      max_length=[1.] * 1,
                      min_deltav_ks91=[-2.5] * 1,
                      max_airmass=[2] * 1,
                      min_twilight_ang=[8] * 1,
                      min_moon_sep=[15] * 1)

    clist.add_cadence(name='single_4x1',
                      nepochs=4,
                      skybrightness=[1.] * 4,
                      delta=[-1.] * 4,
                      delta_min=[-1.] * 4,
                      delta_max=[-1.] * 4,
                      nexp=[1] * 4,
                      max_length=[1.] * 4,
                      min_deltav_ks91=[-2.5] * 4,
                      max_airmass=[2] * 4,
                      min_twilight_ang=[8] * 4,
                      min_moon_sep=[15] * 4)

    clist.add_cadence(name='single_10x1',
                      nepochs=10,
                      skybrightness=[1.] * 10,
                      delta=[-1.] * 10,
                      delta_min=[-1.] * 10,
                      delta_max=[-1.] * 10,
                      nexp=[1] * 10,
                      max_length=[1.] * 10,
                      min_deltav_ks91=[-2.5] * 10,
                      max_airmass=[2] * 10,
                      min_twilight_ang=[8] * 10,
                      min_moon_sep=[15] * 10)

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
    clist.reset()

    clist.add_cadence(name='mwm_tess_rgb_2x1',
                      nepochs=2,
                      skybrightness=[1., 1.],
                      delta=[0., 26.],
                      delta_min=[0., 1.],
                      delta_max=[0., 3000.],
                      nexp=[1, 1],
                      max_length=[1., 1.],
                      min_deltav_ks91=[-2.5, -2.5],
                      max_airmass=[2, 2],
                      min_twilight_ang=[8, 8],
                      min_moon_sep=[15, 15])

    clist.add_cadence(name='csc_faint_boss_1x4',
                      nepochs=1,
                      skybrightness=[0.35],
                      delta=[-1.],
                      delta_min=[-1.],
                      delta_max=[-1.],
                      nexp=[4],
                      max_length=[1.],
                      min_deltav_ks91=[-1.5],
                      max_airmass=[1.4],
                      min_twilight_ang=[15],
                      min_moon_sep=[35])

    ok, epochs_list = clist.cadence_consistency('mwm_tess_rgb_2x1', 'csc_faint_boss_1x4')
    assert ok is False


def test_cadence_consistency_4():
    clist = cadence.CadenceList()
    clist.reset()

    clist.add_cadence(name='dark_2x2',
                      nepochs=2,
                      skybrightness=[0.35, 0.35],
                      delta=[0., 26.],
                      delta_min=[0., 1.],
                      delta_max=[0., 3000.],
                      nexp=[2, 2],
                      max_length=[1., 1.],
                      min_deltav_ks91=[-1.5] * 2,
                      max_airmass=[1.4] * 2,
                      min_twilight_ang=[15] * 2,
                      min_moon_sep=[35] * 2)

    clist.add_cadence(name='dark_6x3',
                      nepochs=6,
                      skybrightness=[0.35] * 6,
                      delta=[0., 26., 26., 26., 26., 26.],
                      delta_min=[0., 1., 1., 1., 1., 1.],
                      delta_max=[0., 3000., 3000., 3000., 3000., 3000.],
                      nexp=[3] * 6,
                      max_length=[1.] * 6,
                      min_deltav_ks91=[-1.5] * 6,
                      max_airmass=[1.4] * 6,
                      min_twilight_ang=[15] * 6,
                      min_moon_sep=[35] * 6)

    clist.add_cadence(name='dark_4x1',
                      nepochs=4,
                      skybrightness=[0.35] * 4,
                      delta=[-1.] * 4,
                      delta_min=[-1.] * 4,
                      delta_max=[-1.] * 4,
                      nexp=[1] * 4,
                      max_length=[1.] * 4,
                      min_deltav_ks91=[-1.5] * 4,
                      max_airmass=[1.4] * 4,
                      min_twilight_ang=[15] * 4,
                      min_moon_sep=[35] * 4)

    clist.add_cadence(name='dark_8x1',
                      nepochs=8,
                      skybrightness=[0.35] * 8,
                      delta=[-1.] * 8,
                      delta_min=[-1.] * 8,
                      delta_max=[-1.] * 8,
                      nexp=[1] * 8,
                      max_length=[1.] * 8,
                      min_deltav_ks91=[-1.5] * 8,
                      max_airmass=[1.4] * 8,
                      min_twilight_ang=[15] * 8,
                      min_moon_sep=[35] * 8)

    ok, epochs_list = clist.cadence_consistency('dark_4x1',
                                                'dark_2x2')
    assert ok is True

    ok, epochs_list = clist.cadence_consistency('dark_2x2',
                                                'dark_6x3')
    assert ok is True

    ok, epochs_list = clist.cadence_consistency('dark_4x1',
                                                'dark_6x3')
    assert ok is True

    ok, epochs_list = clist.cadence_consistency('dark_2x2',
                                                'dark_4x1')
    assert ok is False

    ok, epochs_list = clist.cadence_consistency('dark_8x1',
                                                'dark_2x2')
    assert ok is False

    ok, epochs_list = clist.cadence_consistency('dark_8x1',
                                                'dark_6x3')
    assert ok is True


def test_specific_cadence_consistency_1():
    clist = cadence.CadenceList()
    clist.reset()

    clist.add_cadence(name='dark_2x2',
                      nepochs=2,
                      skybrightness=[0.35, 0.35],
                      delta=[0., 26.],
                      delta_min=[0., 1.],
                      delta_max=[0., 3000.],
                      nexp=[2, 2],
                      max_length=[1., 1.],
                      min_deltav_ks91=[-1.5] * 2,
                      max_airmass=[1.4] * 2,
                      min_twilight_ang=[15] * 2,
                      min_moon_sep=[35] * 2)

    clist.add_cadence(name='dark_6x3',
                      nepochs=6,
                      skybrightness=[0.35] * 6,
                      delta=[0., 26., 26., 26., 26., 26.],
                      delta_min=[0., 1., 1., 1., 1., 1.],
                      delta_max=[0., 3000., 3000., 3000., 3000., 3000.],
                      nexp=[3] * 6,
                      max_length=[1.] * 6,
                      min_deltav_ks91=[-1.5] * 6,
                      max_airmass=[1.4] * 6,
                      min_twilight_ang=[15] * 6,
                      min_moon_sep=[35] * 6)

    clist.add_cadence(name='dark_4x1',
                      nepochs=4,
                      skybrightness=[0.35] * 4,
                      delta=[-1.] * 4,
                      delta_min=[-1.] * 4,
                      delta_max=[-1.] * 4,
                      nexp=[1] * 4,
                      max_length=[1.] * 4,
                      min_deltav_ks91=[-1.5] * 4,
                      max_airmass=[1.4] * 4,
                      min_twilight_ang=[15] * 4,
                      min_moon_sep=[35] * 4)

    clist.add_cadence(name='dark_8x1',
                      nepochs=8,
                      skybrightness=[0.35] * 8,
                      delta=[-1.] * 8,
                      delta_min=[-1.] * 8,
                      delta_max=[-1.] * 8,
                      nexp=[1] * 8,
                      max_length=[1.] * 8,
                      min_deltav_ks91=[-1.5] * 8,
                      max_airmass=[1.4] * 8,
                      min_twilight_ang=[15] * 8,
                      min_moon_sep=[35] * 8)

    add_cadence_single_nxm(n=4, m=1)
    add_cadence_single_nxm(n=3, m=1)
    add_cadence_mixed2_nxm(n=5, m=1)
    add_cadence_mixed2_nxm(n=8, m=1)
    add_cadence_mixed2_nxm(n=15, m=1)
    add_cadence_mixedreverse2_nxm(n=5, m=1)

    ok, epochs_list = clist.specific_cadence_consistency('single_4x1', 'dark_8x1', [0, 2])
    assert ok is True

    ok, epochs_list = clist.specific_cadence_consistency('mixed2_5x1', 'mixed2_5x1',
                                                         [0, 1, 2])
    assert ok is True
    assert epochs_list[0] == [0, 1, 2]

    ok, epochs_list = clist.specific_cadence_consistency('mixed2_8x1', 'mixed2_5x1',
                                                         [0, 1, 2, 3])
    assert ok is True
    assert epochs_list[0] == [0, 1, 2, 3]

    ok, epochs_list = clist.specific_cadence_consistency('mixed2_8x1', 'mixed2_5x1',
                                                         [0, 2, 3])
    assert ok is True
    assert epochs_list[0] == [0, 2, 3]

    ok, epochs_list = clist.specific_cadence_consistency('mixed2_8x1', 'mixed2_5x1',
                                                         [2, 3])
    assert ok is True
    assert epochs_list[0] == [0, 1]

    ok, epochs_list = clist.specific_cadence_consistency('mixed2_8x1', 'single_3x1',
                                                         [0, 2])
    assert ok is False

    ok, epochs_list = clist.specific_cadence_consistency('single_3x1', 'mixed2_8x1',
                                                         [0, 1])
    assert ok is True
    assert epochs_list[0] == [0, 1]

    ok, epochs_list = clist.specific_cadence_consistency('dark_4x1',
                                                         'dark_2x2', [0, 1, 2, 3])
    assert ok is True
    assert epochs_list == [[0, 0, 1, 1]]

    ok, epochs_list = clist.specific_cadence_consistency('dark_8x1',
                                                         'dark_4x1', [0, 2, 4, 7])
    assert ok is True
    assert epochs_list == [[0, 1, 2, 3]]



def test_cadence_evaluate_next():
    clist = cadence.CadenceList()
    clist.reset()

    clist.add_cadence(name='dark_2x2',
                      nepochs=2,
                      skybrightness=[0.35, 0.35],
                      delta=[0., 26.],
                      delta_min=[0., 1.],
                      delta_max=[0., 40.],
                      nexp=[2, 2],
                      max_length=[1., 1.],
                      min_deltav_ks91=[-1.5] * 2,
                      max_airmass=[1.4] * 2,
                      min_twilight_ang=[15] * 2,
                      min_moon_sep=[35] * 2)

    test_cadence = clist.cadences["dark_2x2"]

    idx = 0
    mjd_past = 0  # this would be returned by epochs_completed
    mjd_next = 59900
    skybrightness_next = 0.2
    observable, priority = \
    test_cadence.evaluate_next(epoch_idx=idx,
                               mjd_past=mjd_past,
                               mjd_next=mjd_next,
                               skybrightness_next=skybrightness_next,
                               moon_dist=50,
                               deltaV=-0.8,
                               airmass=1.2)

    assert observable

    skybrightness_next = 0.5
    observable, priority = \
    test_cadence.evaluate_next(epoch_idx=idx,
                               mjd_past=mjd_past,
                               mjd_next=mjd_next,
                               skybrightness_next=skybrightness_next,
                               moon_dist=50,
                               deltaV=-0.8,
                               airmass=1.2)

    assert not observable

    idx = 1
    mjd_past = 59895  # this would be returned by epochs_completed
    mjd_next = 59900
    skybrightness_next = 0.2
    observable, priority1 = \
    test_cadence.evaluate_next(epoch_idx=idx,
                               mjd_past=mjd_past,
                               mjd_next=mjd_next,
                               skybrightness_next=skybrightness_next,
                               moon_dist=50,
                               deltaV=-0.8,
                               airmass=1.2)

    mjd_next = mjd_past + 25
    observable, priority2 = \
    test_cadence.evaluate_next(epoch_idx=idx,
                               mjd_past=mjd_past,
                               mjd_next=mjd_next,
                               skybrightness_next=skybrightness_next,
                               moon_dist=50,
                               deltaV=-0.8,
                               airmass=1.2)

    assert observable
    assert priority2 > priority1

    # no longer valid because delta_max isn't strictly respected anymore
    # idx = 1
    # mjd_past = 59800  # this would be returned by epochs_completed
    # mjd_next = 59900
    # skybrightness_next = 0.2
    # observable, priority = \
    # test_cadence.evaluate_next(epoch_idx=idx,
    #                            mjd_past=mjd_past,
    #                            mjd_next=mjd_next,
    #                            skybrightness_next=skybrightness_next)

    # assert not observable


def test_cadence_consistency_merge():
    clist = cadence.CadenceList()
    clist.reset()

    clist.add_cadence(name='dark_2x2',
                      nepochs=2,
                      skybrightness=[0.35, 0.35],
                      delta=[-1., -1.],
                      delta_min=[-1., -1.],
                      delta_max=[-1., -1.],
                      nexp=[2, 2],
                      max_length=[1.] * 2,
                      min_deltav_ks91=[-1.5] * 2,
                      max_airmass=[1.4] * 2,
                      min_twilight_ang=[15] * 2,
                      min_moon_sep=[35] * 2)

    clist.add_cadence(name='dark_2x4',
                      nepochs=2,
                      skybrightness=[0.35] * 6,
                      delta=[0., 26., 26., 26., 26., 26.],
                      delta_min=[0., 1., 1., 1., 1., 1.],
                      delta_max=[0., 3000., 3000., 3000., 3000., 3000.],
                      nexp=[4] * 2,
                      max_length=[1.] * 2,
                      min_deltav_ks91=[-1.5] * 2,
                      max_airmass=[1.4] * 2,
                      min_twilight_ang=[15] * 2,
                      min_moon_sep=[35] * 2)

    ok, epochs_list, nexps_list = clist.cadence_consistency('dark_2x2',
                                                            'dark_2x4',
                                                            merge_epochs=True)
    assert ok is True

    assert len(nexps_list[0]) == 1
    assert nexps_list[0][0] == 4

    assert len(nexps_list[1]) == 2
    assert nexps_list[1][0] == 2
    assert nexps_list[1][1] == 2

    assert len(nexps_list[2]) == 1
    assert nexps_list[2][0] == 4
