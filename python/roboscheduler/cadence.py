import sys
import numpy as np
import fitsio
from ortools.constraint_solver import pywrapcp

try:
    import sdssdb.peewee.sdss5db.targetdb as targetdb
    _database = True
except:
    _database = False

# Definition to use when writing to ndarray
fits_type = np.dtype('a40')


# Class to define a singleton
class CadenceSingleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(CadenceSingleton,
                                        cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class Cadence(object):
    """Description of a cadence

    Parameters:
    ----------

    nexposures : np.int32
        number of exposures

    lunation : ndarray of np.float32
        maximum lunation for each exposure

    delta : ndarray of np.float32
        desired offset for each exposure from previous (days)

    delta_min : ndarray of np.float32
        minimum delta to allow (days)

    delta_max : ndarray of np.float32
        maximum delta to allow (days)

    instrument : list of str
            instrument for each exposure

    Attributes:
    ----------

    nexposures : np.int32
        number of exposures

    lunation : ndarray of np.float32
        maximum lunation for each exposure

    delta : ndarray of np.float32
        desired offset for each exposure from previous (days)

    delta_min : ndarray of np.float32
        minimum delta to allow (days)

    delta_max : ndarray of np.float32
        maximum delta to allow (days)

    instrument : list of str
            instrument for each exposure

    nepochs : np.int32
        number of epochs

    epoch_indx : ndarray of np.int32
        index into delta for first exposure at each epoch

    epoch_nexposures : ndarray of np.float32
        number of exposures for each separate epoch

    Methods:
    -------

    Notes:
    -----

    The first delta == 0 in all cadences.

    delta == 0 for any exposure except the first means to "bundle" the
    exposure with the previous exposure into a single epoch. delta_min
    and delta_max for such exposures are automatically set to zero on
    input.

    The last epochs in the list may have delta = -1. These will be
    treated as unconstrained. Only single epoch cadences will be allowed.
"""
    def __init__(self, nexposures=None, lunation=None, delta=None,
                 delta_min=None, delta_max=None, instrument=None):
        self.nexposures = np.int32(nexposures)
        self.lunation = np.zeros(self.nexposures, dtype=np.float32) + lunation
        self.instrument = np.array(instrument)
        self.delta = np.zeros(self.nexposures, dtype=np.float32) + delta
        self.delta_min = (np.zeros(self.nexposures, dtype=np.float32) +
                          delta_min)
        self.delta_max = (np.zeros(self.nexposures, dtype=np.float32) +
                          delta_max)
        izero = np.where(self.delta == 0.)[0]
        self.delta_min[izero] = 0.
        self.delta_max[izero] = 0.
        self._create_epochs()
        iapogee = np.where(self.instrument == 'APOGEE')[0]
        if(len(iapogee) > 0):
            self.requires_apogee = 1
        else:
            self.requires_apogee = 0
        iboss = np.where(self.instrument == 'BOSS')[0]
        if(len(iboss) > 0):
            self.requires_boss = 1
        else:
            self.requires_boss = 0
        self.nexposures_pack = len(np.where(self.delta != -1)[0])
        self.nepochs_pack = len(np.where(self.delta[self.epoch_indx] != -1)[0])
        return

    def __str__(self):
        return(self.exposure_text())

    def exposure_text(self):
        """Display list of exposures as text"""
        out = "nexposures={nexposures}\n".format(nexposures=self.nexposures)
        out = out + " lunation = "
        for i in np.arange(self.nexposures):
            out = out + " {l:.3f}".format(l=self.lunation[i])
        out = out + "\n"
        out = out + " delta = "
        for i in np.arange(self.nexposures):
            out = out + " {d:.3f}".format(d=self.delta[i])
        out = out + "\n"
        out = out + " delta_min = "
        for i in np.arange(self.nexposures):
            out = out + " {s:.3f}".format(s=self.delta_min[i])
        out = out + "\n"
        out = out + " delta_max = "
        for i in np.arange(self.nexposures):
            out = out + " {s:.3f}".format(s=self.delta_max[i])
        out = out + "\n"
        out = out + " instrument = "
        for i in np.arange(self.nexposures):
            out = out + " {s}".format(s=self.instrument[i])
        out = out + "\n"
        return(out)

    def epoch_text(self):
        """Display list of epochs as text"""
        epoch_indx = self.epoch_indx
        out = "nepochs={nepochs}\n".format(nepochs=self.nepochs)
        out = out + " nexposures = "
        for i in np.arange(self.nepochs):
            out = out + " {l}".format(l=self.epoch_nexposures[i])
        out = out + "\n"
        out = out + " lunation = "
        for i in np.arange(self.nepochs):
            out = out + " {l:.3f}".format(l=self.lunation[epoch_indx[i]])
        out = out + "\n"
        out = out + " delta = "
        for i in np.arange(self.nepochs):
            out = out + " {d:.3f}".format(d=self.delta[epoch_indx[i]])
        out = out + "\n"
        out = out + " delta_min = "
        for i in np.arange(self.nepochs):
            out = out + " {s:.3f}".format(s=self.delta_min[epoch_indx[i]])
        out = out + "\n"
        out = out + " delta_max = "
        for i in np.arange(self.nepochs):
            out = out + " {s:.3f}".format(s=self.delta_max[epoch_indx[i]])
        out = out + "\n"
        out = out + " instrument = "
        for i in np.arange(self.nepochs):
            out = out + " {s}".format(s=self.instrument[epoch_indx[i]])
        out = out + "\n"
        return(out)

    def _arrayify(self, quantity=None, dtype=np.float64):
        """Cast quantity as ndarray of numpy.float64"""
        try:
            length = len(quantity)
        except TypeError:
            length = 1
        return np.zeros(length, dtype=dtype) + quantity

    def _create_epochs(self):
        """Define epochs based on exposure list"""
        epoch_indx = [0]
        self.nepochs = 1
        for indx in np.arange(self.nexposures - 1) + 1:
            if(self.delta[indx] != 0.):
                epoch_indx.append(indx)
                self.nepochs = self.nepochs + 1
        self.epoch_indx = np.array(epoch_indx, dtype=np.int32)
        self.epoch_nexposures = np.zeros(self.nepochs, dtype=np.int32)
        for indx in np.arange(self.nepochs - 1):
            self.epoch_nexposures[indx] = (self.epoch_indx[indx + 1] -
                                           self.epoch_indx[indx])
        self.epoch_nexposures[-1] = self.nexposures - self.epoch_indx[-1]
        return

    def evaluate_next(self, mjd_past=None, mjd_next=None,
                      lunation_next=None, check_lunation=True):
        """Evaluate next choice of observation (not well-tested)"""
        nexposures_past = len(mjd_past)
        if(nexposures_past >= self.nexposures):
            return(False)
        ok_lunation = ((lunation_next < self.lunation[nexposures_past]) |
                       (check_lunation is False))
        if(nexposures_past == 0):
            return(ok_lunation)
        delta = mjd_next - mjd_past[nexposures_past - 1]
        dlo = self.delta_min[nexposures_past]
        dhi = self.delta_max[nexposures_past]
        return(ok_lunation & (delta >= dlo) & (delta <= dhi))


class Packing(object):
    """Packing of observations into a cadence

    Parameters:
    ----------

    field_cadence : str
        name of cadence to be packing in to

    Attributes:
    ----------

    epoch_targets : list of ndarray of np.int32
        target indices for each epoch

    exposure_targets : array of target indices for each exposure

    Methods:
    -------

    add_target() : add a new cadence
    check_target() : check if a new target solution is allowed
"""
    def __init__(self, field_cadence=None):
        import roboscheduler.cadence2 as cadence
        self.field_cadence = field_cadence
        self.clist = cadence.CadenceList()
        if(self.field_cadence not in self.clist.cadences):
            print("No cadence {fc}.".format(fc=self.field_cadence))
            return
        self.reset()
        return

    def __str__(self):
        return(self.text())

    def text(self):
        out = "Packing of {fc}\n".format(fc=self.field_cadence)
        out = out + " navailable = "
        for i in np.arange(self.nepochs):
            out = out + " {l}".format(l=self.epoch_nexposures[i])
        out = out + "\n"
        out = out + " nused = "
        for i in np.arange(self.nepochs):
            out = out + " {l}".format(l=self.epoch_nused[i])
        out = out + "\n"
        for i in np.arange(self.nepochs):
            out = out + "Epoch #{i}: ".format(i=i + 1)
            for j in np.arange(self.epoch_nexposures[i]):
                out = out + " {x} ".format(x=self.epoch_targets[i][j])
            out = out + "\n"

        return(out)

    def reset(self):
        c = self.clist.cadences[self.field_cadence]
        self.nepochs = c.nepochs
        self.epoch_targets = [np.zeros(x, dtype=np.int32) - 1
                              for x in c.epoch_nexposures]
        self.epoch_nexposures = c.epoch_nexposures
        self.epoch_nused = np.zeros(self.nepochs, dtype=np.int32)
        self.set_exposures()
        return

    def check_target(self, target_cadence=None):
        ok, solns = self.clist.cadence_consistency(target_cadence,
                                                   self.field_cadence,
                                                   return_solutions=True,
                                                   epoch_level=True)
        if(ok is False):
            return(False, [], [])

        for soln in solns:
            pack_soln = soln[0].copy()
            fill_grid = soln[1]
            navail = self.epoch_nexposures - self.epoch_nused
            nneed = self.clist.cadences[target_cadence].epoch_nexposures.copy()
            soln_ok = True
            for indx in np.arange(len(pack_soln), dtype=np.int32):
                pack_epoch = pack_soln[indx]
                if(navail[pack_epoch] < nneed[indx]):
                    soln_ok = False
                    break
                else:
                    navail[indx] = navail[indx] - nneed[indx]
            if(soln_ok):
                if(fill_grid is not None):
                    nexp1 = nneed[self.clist.cadences[target_cadence].nepochs_pack:]
                    nexp2 = navail
                    ok, fill_epoch_targets = self.clist.fill_grid(fill_grid,
                                                                  nexp1=nexp1,
                                                                  nexp2=nexp2)
                    if(ok):
                        return(True, pack_soln, fill_epoch_targets)
                else:
                    return(True, pack_soln, [])

        return(False, [], [])

    def set_exposures(self):
        self.exposures = np.zeros(0, dtype=np.int32)
        n = 0
        for et in self.epoch_targets:
            self.exposures = np.append(self.exposures, et)
            n = n + len(et)
        return

    def add_target(self, target_id=None, target_cadence=None):
        ok, pack_soln, fill_epoch_targets = self.check_target(target_cadence=target_cadence)
        if(ok is False):
            return(False)

        for indx in np.arange(self.clist.cadences[target_cadence].nepochs_pack):
            epoch = pack_soln[indx]
            n = self.epoch_nused[epoch]
            nneed = self.clist.cadences[target_cadence].epoch_nexposures[indx]
            self.epoch_targets[epoch][n:n + nneed] = target_id
            self.epoch_nused[epoch] = n + nneed

        for epoch in np.arange(len(fill_epoch_targets)):
            n = self.epoch_nused[epoch]
            nfill = len(fill_epoch_targets[epoch])
            self.epoch_targets[epoch][n:n + nfill] = target_id
            self.epoch_nused[epoch] = n + nfill

        self.set_exposures()

        return(True)

    def pack_targets_greedy(self, target_ids=None, target_cadences=None,
                            value=None):
        """Pack targets into a cadence greedily.

        Parameters:
        ----------

        target_ids : ndarray of np.int32
            target ID numbers (default np.arange(len(target_cadences)))

        target_cadences : list of strings
            names of the target cadences

        value : ndarray of np.float32
            value for each target (default all 1s)

        Notes:
        -----

        Starts with highest "value" targets, and assigns in a greedy fashion.
        Will not observe partial cadences.
"""
        ntargets = len(target_cadences)

        if(value is None):
            value = np.ones(ntargets)
        else:
            value = np.array(value)

        if(target_ids is None):
            target_ids = np.arange(ntargets, dtype=np.int32)
        else:
            target_ids = np.array(target_ids)

        isort = np.flip(np.argsort(value), axis=0)
        for i in isort:
            self.add_target(target_id=target_ids[i],
                            target_cadence=target_cadences[i])

        self.set_exposures()

        return


class CadenceList(object, metaclass=CadenceSingleton):
    """List of cadences available (singleton)

    Parameters:
    ----------

    Attributes:
    ----------

    ncadences : np.int32, int
         number of different cadences

    cadences : dictionary
         dictionary of Cadence objects

    Methods:
    -------

    reset() : remove all the current cadences
    add_cadence() : add a new cadence
    check_exposures() : are two exposure sets consistent?
    cadence_consistency(): is cadence #1 consistent with cadence #2?
    pack_targets(): pack targets into a cadence optimally
    pack_targets_single(): pack single-shot cadence targets
    pack_targets_greedy(): pack targets into a cadence in a greedy way
    fromarray(): add to cadence list from an ndarray
    fromfits(): add to cadence list from a FITS file
    toarray(): return an ndarray with cadence list
    epoch_array(): return an ndarray with epoch-oriented list of cadences
    todb(): insert cadences into the targetdb
    fromdb(): extract cadences into the targetdb

    Notes:
    -----

    This is a singleton, so there can only be one CadenceList defined
    within any session.

    The cadence_consistency() method stores results for any set of its
    inputs in a dictionary for lookup later; note that if you alter or
    replace a cadence without using reset() then the cached version of
    the cadence consistencies may not be accurate.
"""
    def __init__(self):
        self.reset()
        self.max_nsolns = 100
        return

    def reset(self):
        """Reset cadence list to be empty"""
        self.ncadences = 0
        self.cadences = dict()
        self._cadence_consistency = dict()
        return

    def add_cadence(self, name=None, *args, **kwargs):
        """Add a cadence to the list of cadences

        Parameters:
        ----------

        name : string
            dictionary name of cadence

        nexposures : np.int32
            number of exposures (default 1)

        lunation : ndarray of np.float32
            maximum lunation for each exposure (default [1.])

        delta : ndarray of np.float32
            day for exposure (default [0.])

        delta_min : ndarray of np.float32
            allowance for variation from cadence (default [1.])

        delta_max : ndarray of np.float32
            allowance for variation from cadence (default [1.])

        instrument : list of str
            instrument for each exposure

        Notes:
        -----

        The last epochs in the list may have delta = -1. These will be
        treated as unconstrained. Only single epoch cadences will be allowed.
"""
        cadence = Cadence(*args, **kwargs)
        self.cadences[name] = cadence
        self.ncadences = len(self.cadences.keys())

    def check_exposures(self, one=None, two=None, indx2=None, sub1=None,
                        epoch_level=True, details=False):
        """Is exposure set in cadence two consistent with cadence one?

        Parameters:
        ----------

        one : string
            name of cadence #1

        two : string
            name of cadence #2

        indx2 : ndarray of np.int32
            exposures in cadence #2
        
        sub1 : ndarray of np.int32
            sequence to compare in cadence #1

        epoch_level : boolean
            compare sequences at epoch level not exposure level (default True)

        details : boolean
            if True, output details on failure (default False)

        Returns:
        -------

        ok : boolean
            exposure sequences are compatible

        details : str
            details on result if failure

        Notes:
        -----

        If cadence one contains epochs with delta = -1 they will be ignored.
"""
        eps = 1.e-3  # generously deal with round-off
        nexp = len(indx2)

        if(sub1 is None):
            if(epoch_level):
                sub1 = np.where(self.cadences[one].epoch_delta != -1)[0]
            else:
                sub1 = np.where(self.cadences[one].delta != -1)[0]

        # Check number of exposures, if at epoch level
        if(epoch_level):
            for indx in np.arange(nexp):
                if(self.cadences[one].epoch_nexposures[sub1[indx]] >
                   self.cadences[two].epoch_nexposures[indx2[indx]]):
                    if(details):
                        return(False, "not_enough_exposures")
                    else:
                        return(False)

        # For the subsequent checks, convert to exposure index if we
        # are at the epoch level
        if(epoch_level):
            indx2 = self.cadences[two].epoch_indx[indx2]
            sub1 = self.cadences[one].epoch_indx[sub1]

        # Check lunations
        for indx in np.arange(nexp):
            if(self.cadences[one].lunation[sub1[indx]] <
               self.cadences[two].lunation[indx2[indx]] - eps):
                if(details):
                    return(False, "lunation")
                else:
                    return(False)

        # Check deltas
        for indx in np.arange(nexp - 1) + 1:
            delta1 = self.cadences[one].delta[sub1[indx]]
            dlo1 = self.cadences[one].delta_min[sub1[indx]]
            dhi1 = self.cadences[one].delta_max[sub1[indx]]
            delta2 = self.cadences[two].delta[indx2[indx - 1] + 1:
                                              indx2[indx] + 1].sum()
            dlo2 = self.cadences[two].delta_min[indx2[indx - 1] + 1:
                                                indx2[indx] + 1].sum()
            dhi2 = self.cadences[two].delta_max[indx2[indx - 1] + 1:
                                                indx2[indx] + 1].sum()
            if(delta1 > 0.):  # normal case
                if(dlo1 >= dlo2 + eps):
                    if(details):
                        return(False, "delta_too_soon")
                    else:
                        return(False)
                if(dhi1 <= dhi2 - eps):
                    if(details):
                        return(False, "delta_too_late")
                    else:
                        return(False)
            elif(delta1 == 0.):  # adjacent exposures
                if(indx2[indx] > indx2[indx - 1] + 1):  # must be adjacent
                    if(details):
                        return(False, "exposures_not_adjacent")
                    else:
                        return(False)
                if(delta2 > 0.):  # must be adjacent
                    if(details):
                        return(False, "exposures_not_adjacent")
                    else:
                        return(False)

        if(details):
            return(True, "")
        else:
            return(True)

    def cadence_consistency(self, one, two, return_solutions=True,
                            epoch_level=True):
        """Is cadence #1 consistent with cadence #2?

        Parameters:
        ----------

        one : string
            name of cadence #1

        two : string
            name of cadence #2

        return_solutions: boolean
            return list of solutions? (default False)

        epoch_level : boolean
            compare sequences at epoch level not exposure level (default True)

        Returns:
        -------

        ok : int
            1 if there is a solution, 0 otherwise

        solutions : list (if return_solutions is True)
            list of solutions, see below.

        Notes:
        -----

        Each cadence may have a set of epochs with delta = -1.
        There are nexposures_pack_1 epochs with delta >= 0 in one
        There are nexposures_fill_1 epochs with delta == -1 in one
        There are nexposures_pack_2 epochs with delta >= 0 in two
        There are nexposures_fill_2 epochs with delta == -1 in two

        A solution is a tuple, with
           solns - an array of nexposures_pack_1 exposures within cadence 2
           fill - a 2-d array (nexposures_fill_1, nexposures_2) expressing
                  which exposures can be filled into for this solution
"""
        # Return cached results
        cache_key = (one, two, epoch_level, return_solutions)
        if(cache_key in self._cadence_consistency):
            return(self._cadence_consistency[cache_key])

        if(epoch_level):
            npack1 = self.cadences[one].nepochs_pack
            npack2 = self.cadences[two].nepochs_pack
            nfill1 = self.cadences[one].nepochs - npack1
            nfill2 = self.cadences[two].nepochs - npack2
        else:
            npack1 = self.cadences[one].nexposures_pack
            npack2 = self.cadences[two].nexposures_pack
            nfill1 = self.cadences[one].nexposures - npack1
            nfill2 = self.cadences[two].nexposures - npack2

        possibles = []

        # Special case if same cadence
        if(one == two):
            success = True
            possibles = list(np.arange(npack1))
            ifill1 = np.arange(nfill1)
            ifill2 = np.arange(nfill2)
            fill = np.zeros((nfill1, npack2 + nfill2), dtype=np.bool)
            for i1 in ifill1:
                for i2 in ifill2:
                    fill[i1, npack2 + i2] = self.check_exposures(one=one, two=two,
                                                                 indx2=[npack2 + i2],
                                                                 sub1=[npack1 + i1],
                                                                 epoch_level=epoch_level)
            if(return_solutions):
                self._cadence_consistency[cache_key] = (success,
                                                        [(possibles, fill)])
            else:
                self._cadence_consistency[cache_key] = success
            return(self._cadence_consistency[cache_key])

        if(npack1 == 0):
            possibles = [[]]
        else:
            # Check which exposures you can start on
            for first in np.arange(npack2 - npack1 + 1):
                ok = self.check_exposures(one=one, two=two, indx2=[first],
                                          sub1=[0], epoch_level=epoch_level)
                if(ok):
                    possibles.append([first])
                if(len(possibles) == 0):
                    success = False
                    if(return_solutions):
                        self._cadence_consistency[cache_key] = (success,
                                                                possibles)
                    else:
                        self._cadence_consistency[cache_key] = success
                    return(self._cadence_consistency[cache_key])

            # Now find sequences starting from there
            for nsub1 in np.arange(npack1 - 1) + 2:
                current_possibles = possibles
                possibles = []
                for indx in range(len(current_possibles)):
                    possible = current_possibles[indx]
                    remaining_start = possible[-1] + 1
                    nremaining = npack2 - possible[-1] - 1
                    ok = 1
                    if(nremaining >= npack1 - len(possible)):
                        for next_possible in (remaining_start +
                                              np.arange(nremaining)):
                            try_possible = possible.copy()
                            try_possible.append(next_possible)
                            ok = self.check_exposures(one=one, two=two,
                                                      indx2=try_possible,
                                                      sub1=np.arange(nsub1),
                                                      epoch_level=epoch_level)
                            if(ok):
                                possibles.append(try_possible)

                if(len(possibles) == 0):
                    success = False
                    if(return_solutions):
                        self._cadence_consistency[cache_key] = (success,
                                                                possibles)
                    else:
                        self._cadence_consistency[cache_key] = success
                    return(self._cadence_consistency[cache_key])

        if((len(possibles) > 0) & (nfill1 > 0)):
            pack_possibles = possibles.copy()
            possibles = []
            for pack_possible in pack_possibles:
                fill = np.zeros((nfill1, npack2 + nfill2), dtype=np.bool)
                ifill1 = np.arange(nfill1)
                ifill2 = np.arange(npack2 + nfill2)
                for i1 in ifill1:
                    for i2 in np.arange(npack2 + nfill2, dtype=np.int32):
                        if(i2 not in pack_possible):
                            fill[i1, i2] = self.check_exposures(one=one,
                                                                two=two,
                                                                indx2=[i2],
                                                                sub1=[i1 + npack1],
                                                                epoch_level=epoch_level)
                if(epoch_level):
                    nexp1 = self.cadences[one].epoch_nexposures[ifill1]
                    nexp2 = self.cadences[two].epoch_nexposures[ifill2]
                else:
                    nexp1 = self.cadences[one].nexposures[ifill1]
                    nexp2 = self.cadences[two].nexposures[ifill2]

                ok = self.fill_grid(fill=fill, nexp1=nexp1, nexp2=nexp2)
                if(ok):
                    possibles.append((pack_possible, fill))
        else:
            pack_possibles = possibles.copy()
            possibles = []
            for pack_possible in pack_possibles:
                possibles.append((pack_possible, None))

        success = len(possibles) > 0
        if(return_solutions):
            self._cadence_consistency[cache_key] = (success, possibles)
        else:
            self._cadence_consistency[cache_key] = success
        return(self._cadence_consistency[cache_key])

    def fill_grid(self, fill=None, nexp1=None, nexp2=None):
        """Checks if a grid can be filled

        Parameters:
        ----------

        fill : 2-D ndarray of bool
            whether epoch of second dimension can satisfy epoch of first

        nexp1 : ndarray of np.int32
            number of exposures in first dimension

        nexp2 : ndarray of np.int32
            number of exposures in second dimension

        Returns:
        -------

        ok : bool
            True if the grid can be filled, False if not
"""
        n1 = fill.shape[0]
        n2 = fill.shape[1]

        solver = pywrapcp.Solver("fill_grid")

        gridvars = []
        for indx1 in np.arange(n1):
            indx1vars = []
            for indx2 in np.arange(n2):
                nfill = int(0)
                if(fill[indx1, indx2]):
                    nfill = int(nexp2[indx2])
                name = "{indx1}-{indx2}".format(indx1=indx1, indx2=indx2)
                tmpvar = solver.IntVar(0, nfill, name)
                indx1vars.append(tmpvar)
            gridvars.append(indx1vars)

        for indx1 in np.arange(n1):
            indx1vars = gridvars[indx1]
            solver.Add(solver.Sum(indx1vars) == int(nexp1[indx1]))

        for indx2 in np.arange(n2):
            indx2vars = [x[indx2] for x in gridvars]
            solver.Add(solver.Sum(indx2vars) <= int(nexp2[indx2]))

        allvars = [var for indx1vars in gridvars
                   for var in indx1vars]

        db = solver.Phase(allvars, solver.CHOOSE_FIRST_UNBOUND,
                          solver.ASSIGN_MIN_VALUE)

        # Create a solution collector.
        collector = solver.FirstSolutionCollector()

        # Add the decision variables.
        for allvar in allvars:
            collector.Add(allvar)

        tl = solver.TimeLimit(100)
        status = solver.Solve(db, [collector, tl])
        if(status is False):
            return(False, [])

        # Retrieve list of targets for each epoch
        if collector.SolutionCount() > 0:
            epoch_targets = [np.zeros(0, dtype=np.int32)] * n2
            for indx1 in np.arange(n1):
                indx1vars = gridvars[indx1]
                for i2 in range(len(indx1vars)):
                    var = indx1vars[i2]
                    if(collector.Value(0, var)):
                        for i in np.arange(nexp1[indx1]):
                            epoch_targets[i2] = np.append(epoch_targets[i2],
                                                          np.int32(indx1))
            return(True, epoch_targets)
        else:
            return(False, [])

    def fromarray(self, cadences_array=None, nathan=False):
        """Add cadences to ccadence list from an array

        Parameters:
        -----------

        cadences_array : ndarray
            ndarray with columns 'NEXPOSURES', 'LUNATION', 'DELTA',
            'DELTA_MIN', 'DELTA_MAX', 'CADENCE', 'INSTRUMENT'

        nathan : bool
            False if normal format, True if Nathan De Lee format
"""
        col = {'NEXPOSURES': 'NEXPOSURES',
               'LUNATION': 'LUNATION',
               'DELTA': 'DELTA',
               'DELTA_MIN': 'DELTA_MIN',
               'DELTA_MAX': 'DELTA_MAX',
               'INSTRUMENT': 'INSTRUMENT',
               'CADENCE': 'CADENCE'}
        if(nathan is True):
            for k in col.keys():
                col[k] = col[k].lower()
            col['NEXPOSURES'] = 'nepochs'
            for indx in np.arange(len(cadences_array), dtype=np.int32):
                if(cadences_array[indx][col['NEXPOSURES']] > 1):
                    iz = np.where(cadences_array[indx][col['DELTA']] == 0.)[0]
                    cadences_array[indx][col['DELTA_MIN']][iz] = 0.
                    cadences_array[indx][col['DELTA_MAX']][iz] = 0.

        for ccadence in cadences_array:
            nexp = ccadence[col['NEXPOSURES']]
            if(isinstance(ccadence[col['LUNATION']],
                          type(np.zeros(0, dtype=np.float32)))):
                instruments = np.array([ii.decode().strip()
                                        for ii in ccadence[col['INSTRUMENT']][0:nexp]])
                self.add_cadence(nexposures=ccadence[col['NEXPOSURES']],
                                 lunation=ccadence[col['LUNATION']][0:nexp],
                                 delta=ccadence[col['DELTA']][0:nexp],
                                 delta_min=ccadence[col['DELTA_MIN']][0:nexp],
                                 delta_max=ccadence[col['DELTA_MAX']][0:nexp],
                                 name=ccadence[col['CADENCE']].decode().strip(),
                                 instrument=instruments)
            else:
                instruments = np.array([ccadence[col['INSTRUMENT']].decode().strip()])
                self.add_cadence(nexposures=ccadence[col['NEXPOSURES']],
                                 lunation=ccadence[col['LUNATION']],
                                 delta=ccadence[col['DELTA']],
                                 delta_min=ccadence[col['DELTA_MIN']],
                                 delta_max=ccadence[col['DELTA_MAX']],
                                 name=ccadence[col['CADENCE']].decode().strip(),
                                 instrument=instruments)
        return

    def fromfits(self, filename=None, nathan=False):
        """Add cadences to ccadence list from a FITS file

        Parameters:
        -----------

        filename : str
            File name to read from

        Notes:
        -----

        Expects a valid FITS file with columns 'NEXPOSURES',
            'LUNATION', 'DELTA', 'DELTA_MIN', 'DELTA_MAX', 'CADENCE',
            'INSTRUMENT'
"""
        self.cadences_fits = fitsio.read(filename)
        self.fromarray(self.cadences_fits, nathan=nathan)
        return

    def toarray(self):
        """Return cadences as a record array

        Returns:
        -------

        cadences : ndarray
            information on each cadence
"""
        nexps = np.array([c.nexposures for c in self.cadences.values()])
        max_nexp = nexps.max()
        cadence0 = [('CADENCE', fits_type),
                    ('NEXPOSURES', np.int32),
                    ('DELTA', np.float64, max_nexp),
                    ('LUNATION', np.float32, max_nexp),
                    ('DELTA_MAX', np.float32, max_nexp),
                    ('DELTA_MIN', np.float32, max_nexp),
                    ('INSTRUMENT', np.dtype('a10'), max_nexp)]
        cads = np.zeros(self.ncadences, dtype=cadence0)
        names = self.cadences.keys()
        for indx, name in zip(np.arange(self.ncadences), names):
            nexp = self.cadences[name].nexposures
            cads['CADENCE'][indx] = name
            cads['NEXPOSURES'][indx] = nexp
            cads['DELTA'][indx, 0:nexp] = self.cadences[name].delta
            cads['DELTA_MIN'][indx, 0:nexp] = self.cadences[name].delta_min
            cads['DELTA_MAX'][indx, 0:nexp] = self.cadences[name].delta_max
            cads['LUNATION'][indx, 0:nexp] = self.cadences[name].lunation
            cads['INSTRUMENT'][indx, 0:nexp] = self.cadences[name].instrument
        return(cads)

    def epoch_array(self):
        """Return cadence epochs as a record array

        Returns:
        -------

        cadences : ndarray
            information on each cadence
"""
        neps = np.array([c.nepochs for c in self.cadences.values()])
        max_nep = neps.max()
        cadence0 = [('CADENCE', fits_type),
                    ('NEPOCHS', np.int32),
                    ('NEXPOSURES', np.int32, max_nep),
                    ('DELTA', np.float64, max_nep),
                    ('LUNATION', np.float32, max_nep),
                    ('DELTA_MAX', np.float32, max_nep),
                    ('DELTA_MIN', np.float32, max_nep),
                    ('INSTRUMENT', np.dtype('a10'), max_nep)]
        cads = np.zeros(self.ncadences, dtype=cadence0)
        names = self.cadences.keys()
        for indx, name in zip(np.arange(self.ncadences), names):
            nep = self.cadences[name].nepochs
            cads['CADENCE'][indx] = name
            cads['NEPOCHS'][indx] = nep
            epoch_indx = self.cadences[name].epoch_indx
            cads['NEXPOSURES'][indx, 0:nep] = self.cadences[name].epoch_nexposures
            cads['DELTA'][indx, 0:nep] = self.cadences[name].delta[epoch_indx]
            cads['DELTA_MIN'][indx, 0:nep] = self.cadences[name].delta_min[epoch_indx]
            cads['DELTA_MAX'][indx, 0:nep] = self.cadences[name].delta_max[epoch_indx]
            cads['LUNATION'][indx, 0:nep] = self.cadences[name].lunation[epoch_indx]
            instruments = [self.cadences[name].instrument[i] for i in epoch_indx]
            cads['INSTRUMENT'][indx][0:nep] = instruments
        return(cads)

    def todb(self):
        """Insert all cadences into the targetdb
"""
        if(_database is False):
            print("No database available.")
            return()

        # Create dictionary to look up spectrograph pk from instrument name
        spectrographs = targetdb.Spectrograph.select().dicts()
        spectrograph_pk = dict()
        for spectrograph in spectrographs:
            spectrograph_pk[spectrograph['label']] = spectrograph['pk']

        pkdict = targetdb.TargetCadence.select(targetdb.TargetCadence.pk).dicts()
        pks = np.array([p['pk'] for p in pkdict])
        newpks = pks.max() + 1 + np.arange(len(self.cadences))

        for cadence, pk in zip(self.cadences, newpks):
            nexposures = int(self.cadences[cadence].nexposures)
            delta = [float(n) for n in self.cadences[cadence].delta]
            delta_min = [float(n) for n in self.cadences[cadence].delta_min]
            delta_max = [float(n) for n in self.cadences[cadence].delta_max]
            lunation = [float(n) for n in self.cadences[cadence].lunation]
            spectrograph = [spectrograph_pk[n]
                            for n in self.cadences[cadence].instrument]
            targetdb.TargetCadence.insert(pk=pk, name=cadence,
                                          nexposures=nexposures,
                                          delta=delta, lunation=lunation,
                                          delta_min=delta_min,
                                          delta_max=delta_max,
                                          spectrograph_pk=spectrograph).execute()

    def updatedb(self):
        """Update the cadences in the targetdb by name
"""
        if(_database is False):
            print("No database available.")
            return()

        # Create dictionary to look up spectrograph pk from instrument name
        spectrographs = targetdb.Spectrograph.select().dicts()
        spectrograph_pk = dict()
        for spectrograph in spectrographs:
            spectrograph_pk[spectrograph['label']] = spectrograph['pk']

        for cadence in self.cadences:
            update_dict = {targetdb.TargetCadence.nexposures:
                           self.cadences[cadence].nexposures,
                           targetdb.TargetCadence.delta:
                           [float(n) for n in self.cadences[cadence].delta],
                           targetdb.TargetCadence.delta_min:
                           [float(n)
                            for n in self.cadences[cadence].delta_min],
                           targetdb.TargetCadence.delta_max:
                           [float(n)
                            for n in self.cadences[cadence].delta_max],
                           targetdb.TargetCadence.lunation:
                           [float(n) for n in self.cadences[cadence].lunation],
                           targetdb.TargetCadence.spectrograph_pk:
                           [spectrograph_pk[n]
                            for n in self.cadences[cadence].instrument]}
            targetdb.TargetCadence.update(update_dict). \
                where(targetdb.TargetCadence.name == cadence).execute()
        return

    def fromdb(self):
        """Extract cadences into the targetdb
"""
        if(_database is False):
            print("No database available.")
            return()

        # Create dictionary to look up spectrograph pk from instrument name
        spectrographs = targetdb.Spectrograph.select().dicts()
        instrument = dict()
        for spectrograph in spectrographs:
            instrument[spectrograph['pk']] = spectrograph['label']

        cadences = targetdb.TargetCadence.select().dicts()
        for cadence in cadences:
            instruments = np.array([instrument[pk]
                                    for pk in cadence['spectrograph_pk']])
            if(len(instruments) == 0):
                print("No instruments, defaulting to APOGEE")
                instruments = ['APOGEE'] * cadence['nexposures']
            self.add_cadence(name=cadence['name'],
                             nexposures=cadence['nexposures'],
                             delta=np.array(cadence['delta']),
                             delta_min=np.array(cadence['delta_min']),
                             delta_max=np.array(cadence['delta_max']),
                             lunation=np.array(cadence['lunation']),
                             instrument=instruments)
