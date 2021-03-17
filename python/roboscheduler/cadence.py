import pickle
import numpy as np
import fitsio
import roboscheduler.cCadenceCore as cCadenceCore

Instrument = cCadenceCore.Instrument

try:
    import sdssdb.peewee.sdss5db.targetdb as targetdb
    _database = True
except:
    _database = False


def basename(cadence):
    return("_".join(cadence.split('-')[0].split('_')[0:-1]))


def _instrument_name(instrument):
    instrument_name = 'UNKNOWN'
    if(instrument == Instrument.BossInstrument):
        instrument_name = 'BOSS'
    if(instrument == Instrument.ApogeeInstrument):
        instrument_name = 'APOGEE'
    return(instrument_name)


def _name_instrument(name):
    name_instrument = None
    if(name == 'BOSS'):
        name_instrument = Instrument.BossInstrument
    if(name == 'APOGEE'):
        name_instrument = Instrument.ApogeeInstrument
    return(name_instrument)


# Class to define a singleton
class CadenceListSingleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(CadenceListSingleton,
                                        cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class Cadence(cCadenceCore.CadenceCore):
    """Description of a cadence

    Parameters:
    ----------

    name : str
        name of cadence

    nepochs : np.int32
        number of epochs

    instrument : str
            instrument

    skybrightness : ndarray of np.float32
        maximum sky brightness for each exposure

    delta : ndarray of np.float32
        desired offset for each exposure from previous (days)

    delta_min : ndarray of np.float32
        minimum delta to allow (days)

    delta_max : ndarray of np.float32
        maximum delta to allow (days)

    nexp : ndarray of np.int32
        number of exposures each epoch

    Attributes:
    ----------

    nexposures : np.int32
        number of exposures

    skybrightness : ndarray of np.float32
        maximum sky brightness for each exposure

    delta : ndarray of np.float32
        desired offset for each exposure from previous (days)

    delta_min : ndarray of np.float32
        minimum delta to allow (days)

    delta_max : ndarray of np.float32
        maximum delta to allow (days)

    instrument : str
        instrument ('APOGEE' or 'BOSS')

    nepochs : np.int32
        number of epochs

    epoch_indx : ndarray of np.int32
        index into delta for first exposure at each epoch

    epoch_nexposures : ndarray of np.int32
        number of exposures for each separate epoch

    epochs : ndarray of np.int32
        epoch number for each exposure

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
    def __init__(self, name=None, nepochs=None, skybrightness=None,
                 delta=None, delta_min=None, delta_max=None, nexp=None,
                 instrument=None):
        nepochs = np.int32(nepochs)
        skybrightness = self._arrayify(skybrightness, dtype=np.float32)
        if(instrument == 'APOGEE'):
            instrument = cCadenceCore.ApogeeInstrument
        else:
            instrument = cCadenceCore.BossInstrument
        delta = self._arrayify(delta, dtype=np.float32)
        delta_min = self._arrayify(delta_min, dtype=np.float32)
        delta_max = self._arrayify(delta_max, dtype=np.float32)
        nexp = self._arrayify(nexp, dtype=np.int32)
        epoch_indx = np.zeros(nepochs + 1, dtype=np.int32)
        epochs = np.zeros(nexp.sum(), dtype=np.int32)
        super().__init__(name, nepochs, instrument,
                         skybrightness, delta, delta_min,
                         delta_max, nexp, epoch_indx, epochs)
        return

    def as_cadencecore(self):
        cc = cCadenceCore.CadenceCore(self.name, self.nepochs, self.instrument,
                                      self.skybrightness, self.delta,
                                      self.delta_min, self.delta_max, self.nexp,
                                      self.epoch_indx, self.epochs)
        return(cc)

    def _arrayify(self, quantity=None, dtype=np.float64):
        """Cast quantity as ndarray of numpy.float64"""
        try:
            length = len(quantity)
        except TypeError:
            length = 1
        return np.zeros(length, dtype=dtype) + quantity

    def smart_epoch_nexp(self, mjd_past, tolerance=45):
        """Calculate # of observed epochs, allowing
        for more exposures than planned.

        tolerance is in minutes
        """
        nexposures_past = len(mjd_past)
        if(nexposures_past >= self.nexposures):
            return 1

        tolerance = tolerance / 60. / 24.
        obs_epochs = 0
        prev = 0
        for m in mjd_past:
            delta = m - prev
            prev = m
            if delta < tolerance:
                continue
            else:
                obs_epochs += 1

        if obs_epochs >= self.nepochs:
            assert nexposures_past >= self.nexposures, "skipped some exposures!!"
            return 1

        nexposures_next = self.epoch_nexposures[obs_epochs]

        return nexposures_next

    def next_epoch_nexp(self, mjd_past):
        """get number of exposures for elligible epoch"""
        nexposures_past = len(mjd_past)
        if(nexposures_past >= self.nexposures):
            return 1
        epoch_indx = np.where(self.epoch_indx == nexposures_past)[0]
        nexposures_next = self.epoch_nexposures[epoch_indx]
        assert len(nexposures_next) == 1, "epoch selection failed"

        return nexposures_next

    def skybrightness_check(self, mjd_past, skybrightness_next):
        """check lunation for mjd_past against lunation_next"""
        nexposures_past = len(mjd_past)
        if(nexposures_past >= self.nexposures):
            return skybrightness_next <= self.skybrightness[-1]
        return skybrightness_next <= self.skybrightness[nexposures_past]

    def evaluate_next(self, mjd_past=None, mjd_next=None,
                      skybrightness_next=None, check_skybrightness=True,
                      ignoreMax=False):
        """Evaluate next choice of observation

           Returns whether cadence is ok AND how long until
           cadence will become impossible (deltaMax - delta)
        """
        nexposures_past = len(mjd_past)
        if(nexposures_past >= self.nexposures):
            print("done!")
            return(False, 0)

        ok_skybrightness = (self.skybrightness_check(mjd_past, skybrightness_next)|
                       (check_skybrightness is False))
        if(nexposures_past == 0):
            return(ok_skybrightness, 1e6)

        delta = mjd_next - mjd_past[nexposures_past - 1]
        dlo = self.delta_min[nexposures_past]
        dhi = self.delta_max[nexposures_past]
        if(dlo == -1):
            return(ok_skybrightness, 1e6)
        # print("delta {} dhi {} dlo {}".format(delta, dhi, dlo))
        if ignoreMax:
            return(ok_skybrightness & (delta >= dlo), np.abs(dhi - delta))
        return(ok_skybrightness & (delta >= dlo) & (delta <= dhi), dhi - delta)


class CadenceList(object, metaclass=CadenceListSingleton):
    """List of cadences available (singleton)

    Parameters:
    ----------

    Attributes:
    ----------

    ncadences : np.int32, int
         number of different cadences

    cadences : dict
         dictionary of Cadence objects

    Methods:
    -------

    reset() : remove all the current cadences
    add_cadence() : add a new cadence
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

        nepochs : np.int32
            number of epochs (default 1)

        instrument : str
            instrument to use

        skybrightness : ndarray of np.float32
            maximum sky brightness for each exposure (default [1.])

        delta : ndarray of np.float32
            day for exposure (default [0.])

        delta_min : ndarray of np.float32
            allowance for variation from cadence (default [1.])

        delta_max : ndarray of np.float32
            allowance for variation from cadence (default [1.])

        nexp : ndarray of np.int32
            number of exposures per cadence

        Notes:
        -----

        The last epochs in the list may have delta = -1. These will be
        treated as unconstrained. Only single epoch cadences will be allowed.
"""
        cadence = Cadence(name=name, *args, **kwargs)
        self.cadences[name] = cadence
        self.ncadences = len(self.cadences.keys())
        return

    def cadence_consistency(self, one, two, return_solutions=True,
                            epoch_level=True, merge_epochs=False):
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
            [ignores this]

        merge_epochs : bool
            if epochs in the solutions are repeated, merge them

        Returns:
        -------

        ok : int
            1 if there is a solution, 0 otherwise

        solutions : list of ndarrays of np.int32
            list of solutions, where each solution is an array
        
        nexps : list of ndarrays of np.int32
            if epochs are merged, number of exposures for each epoch listed

        Notes:
        -----

        When epochs are merged with merge_epochs you need to keep track
        of how many exposures you are putting into each epoch. In this 
        case 'nexps' is returned.
"""
        cache_key = (one, two, epoch_level, return_solutions, merge_epochs)
        if(cache_key in self._cadence_consistency):
            return(self._cadence_consistency[cache_key])

        onecc = self.cadences[one].as_cadencecore()
        possibles = self.cadences[two].cadence_consistency(onecc)
        success = len(possibles) > 0

        if(return_solutions):
            if(merge_epochs):
                possibles_merge = []
                nexps_merge = []
                for possible in possibles:
                    possible_merge, possible_inverse = np.unique(possible,
                                                                 return_inverse=True)
                    nexp_one = self.cadences[one].nexp
                    nexp_merge = np.zeros(len(possible_merge), dtype=np.int32)
                    for i, nexp in zip(possible_inverse, nexp_one):
                        nexp_merge[i] = nexp_merge[i] + nexp

                    possibles_merge.append(possible_merge)
                    nexps_merge.append(nexp_merge)


                self._cadence_consistency[cache_key] = (success,
                                                        possibles_merge,
                                                        nexps_merge)
            else:
                self._cadence_consistency[cache_key] = (success, possibles)
        else:
            self._cadence_consistency[cache_key] = success

        return(self._cadence_consistency[cache_key])

    def exposure_consistency(self, one, two, iexp):
        """Is cadence #1 consistent with a set of exposures in cadence #2?

        Parameters:
        ----------

        one : string
            name of cadence #1

        two : string
            name of cadence #2

        iexp : ndarray of np.int32
            indices of exposures to check

        Returns:
        -------

        ok : bool
            True if these exposures under cadence #2 satisfy cadence #1
"""
        # Are there any ways cadence one fits into cadence two?
        ok, epochs_list = self.cadence_consistency(one, two)
        if(not ok):
            return False

        cadence_two = self.cadences[two]
        cadence_one = self.cadences[one]

        # Now count how many exposures there are at each epoch for this
        # array of iexp
        epochs_two_got, nexps_got = np.unique(cadence_two.epochs[iexp],
                                              return_counts=True)
        ngot = np.zeros(cadence_two.nepochs, dtype=np.int32)
        ngot[epochs_two_got] = nexps_got

        # Check each solution and see if the current assignment satisfies

        nneed = np.zeros(cadence_two.nepochs, dtype=np.int32)
        for epochs_two in epochs_list:

            # For this solution, count how many exposures are required
            # at each epoch of cadence two
            nneed[:] = 0
            bad = False
            for epoch_one, epoch_two in enumerate(epochs_two):
                nneed[epoch_two] = nneed[epoch_two] + cadence_one.nexp[epoch_one]
                if(nneed[epoch_two] > ngot[epoch_two]):
                    bad = True
                    break

            # For each exposure in the solution, are there enough exposures?
            if(not bad):
                nbad = (ngot < nneed).sum()
                if(nbad == 0):
                    return True

        # If we have gotten here no working solution was found
        return False

    def fromarray(self, cadences_array=None):
        """Add cadences to ccadence list from an array

        Parameters:
        -----------

        cadences_array : ndarray
            ndarray with columns 'NEPOCHS', 'SKYBRIGHTNESS', 'DELTA',
            'DELTA_MIN', 'DELTA_MAX', 'NEXP', 'CADENCE', 'INSTRUMENT'
"""
        for ccadence in cadences_array:
            nepochs = ccadence['NEPOCHS']
            instrument = ccadence['INSTRUMENT']

            self.add_cadence(name=ccadence['CADENCE'],
                             nepochs=ccadence['NEPOCHS'],
                             skybrightness=ccadence['SKYBRIGHTNESS'][0:nepochs],
                             delta=ccadence['DELTA'][0:nepochs],
                             delta_min=ccadence['DELTA_MIN'][0:nepochs],
                             delta_max=ccadence['DELTA_MAX'][0:nepochs],
                             nexp=ccadence['NEXP'][0:nepochs],
                             instrument=instrument)
        return

    def fromfits(self, filename=None, unpickle=False):
        """Add cadences to ccadence list from a FITS file

        Parameters:
        -----------

        filename : str
            File name to read from

        Notes:
        -----

        Expects a valid FITS file with columns 'NEPOCHS',
            'SKYBRIGHTNESS', 'DELTA', 'DELTA_MIN', 'DELTA_MAX', 'CADENCE',
            'NEXP', 'INSTRUMENT'
        """
        self.cadences_fits = fitsio.read(filename)
        self.fromarray(self.cadences_fits)
        if(unpickle):
            fp = open(filename + '.pkl', 'rb')
            cc = pickle.load(fp)
            for key in cc:
                self._cadence_consistency[key] = cc[key]
        return

    def toarray(self):
        """Return cadences as a record array

        Returns:
        -------

        cadences : ndarray
            information on each cadence
        """
        nepochs = np.array([c.nepochs for c in self.cadences.values()])
        max_nexp = nepochs.max()
        cadence0 = [('CADENCE', np.unicode_, 40),
                    ('NEPOCHS', np.int32),
                    ('DELTA', np.float64, max_nexp),
                    ('SKYBRIGHTNESS', np.float32, max_nexp),
                    ('DELTA_MAX', np.float32, max_nexp),
                    ('DELTA_MIN', np.float32, max_nexp),
                    ('NEXP', np.int32, max_nexp),
                    ('INSTRUMENT', np.unicode_, 40)]
        cads = np.zeros(self.ncadences, dtype=cadence0)
        names = self.cadences.keys()
        for indx, name in enumerate(names):
            nepochs = self.cadences[name].nepochs
            cads['CADENCE'][indx] = name
            cads['NEPOCHS'][indx] = nepochs
            cads['DELTA'][indx, 0:nepochs] = self.cadences[name].delta
            cads['DELTA_MIN'][indx, 0:nepochs] = self.cadences[name].delta_min
            cads['DELTA_MAX'][indx, 0:nepochs] = self.cadences[name].delta_max
            cads['NEXP'][indx, 0:nepochs] = self.cadences[name].nexp
            cads['SKYBRIGHTNESS'][indx, 0:nepochs] = self.cadences[name].skybrightness
            cads['INSTRUMENT'][indx] = _instrument_name(self.cadences[name].instrument)
        return(cads)

    def todb(self):
        """Insert all cadences into the targetdb
"""

        if(_database is False):
            print("No database available.")
            return()

        # Create dictionary to look up instrument pk from instrument name
        instruments = targetdb.Instrument.select().dicts()
        instrument_pk = dict()
        for instrument in instruments:
            instrument_pk[instrument['label']] = instrument['pk']

        pkdict = targetdb.Cadence.select(targetdb.Cadence.pk).dicts()
        pks = np.array([p['pk'] for p in pkdict])
        newpks = pks.max() + 1 + np.arange(len(self.cadences))

        for cadence, pk in zip(self.cadences, newpks):
            nexposures = int(self.cadences[cadence].nexposures)
            delta = [float(n) for n in self.cadences[cadence].delta]
            delta_min = [float(n) for n in self.cadences[cadence].delta_min]
            delta_max = [float(n) for n in self.cadences[cadence].delta_max]
            skybrightness = [float(n) for n in self.cadences[cadence].skybrightness]
            instrument = [instrument_pk[n]
                          for n in self.cadences[cadence].instrument]
            targetdb.Cadence.insert(pk=pk, label=cadence,
                                    nexposures=nexposures,
                                    delta=delta, skybrightness=skybrightness,
                                    delta_min=delta_min,
                                    delta_max=delta_max,
                                    instrument_pk=instrument).execute()

    def updatedb(self):
        """Update the cadences in the targetdb by name
"""
        if(_database is False):
            print("No database available.")
            return()

        # Create dictionary to look up instrument pk from instrument name
        instruments = targetdb.Instrument.select().dicts()
        instrument_pk = dict()
        for instrument in instruments:
            instrument_pk[instrument['label']] = instrument['pk']

        for cadence in self.cadences:
            update_dict = {targetdb.Cadence.nexposures:
                           self.cadences[cadence].nexposures,
                           targetdb.Cadence.delta:
                           [float(n) for n in self.cadences[cadence].delta],
                           targetdb.Cadence.delta_min:
                           [float(n)
                            for n in self.cadences[cadence].delta_min],
                           targetdb.Cadence.delta_max:
                           [float(n)
                            for n in self.cadences[cadence].delta_max],
                           targetdb.Cadence.skybrightness:
                           [float(n) for n in self.cadences[cadence].skybrightness],
                           targetdb.Cadence.instrument_pk:
                           [instrument_pk[n]
                            for n in self.cadences[cadence].instrument]}
            targetdb.Cadence.update(update_dict). \
                where(targetdb.Cadence.label == cadence).execute()
        return

    def fromdb(self):
        """Extract cadences into the targetdb
"""
        if(_database is False):
            print("No database available.")
            return()

        # Create dictionary to look up instrument pk from instrument name
        instruments = targetdb.Instrument.select().dicts()
        instrument = dict()
        for cinstrument in instruments:
            instrument[cinstrument['pk']] = cinstrument['label']

        cadences = targetdb.Cadence.select().dicts()

        for cadence in cadences:
            if(cadence['delta'] is None):
                continue

            instruments = np.array([instrument[pk]
                                    for pk in cadence['instrument_pk']])
            if(len(instruments) == 0):
                print("No instruments, defaulting to APOGEE")
                instruments = ['APOGEE'] * cadence['nexposures']
            self.add_cadence(name=cadence['label'],
                             nexposures=cadence['nexposures'],
                             delta=np.array(cadence['delta']),
                             delta_min=np.array(cadence['delta_min']),
                             delta_max=np.array(cadence['delta_max']),
                             skybrightness=np.array(cadence['skybrightness']))
