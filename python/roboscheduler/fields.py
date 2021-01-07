import numpy as np
import fitsio
import roboscheduler.cadence
import sys

# Class to define a singleton
class FieldsSingleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(FieldsSingleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class Fields(object, metaclass=FieldsSingleton):
    """List of fields

    Parameters:
    ----------

    Attributes:
    ----------

    nfields : np.int32
        number of fields

    racen : ndarray of np.float64
        right ascension of each design (J2000 deg)

    deccen : ndarray of np.float64
        declination of each design (J2000 deg)

    fieldid : ndarray of np.int32
        id for each field

    cadence : string
        cadence for each field

    observations : list of ndarray of np.int32
        for each field, indices of observations
"""
    def __init__(self):
        self.nfields = 0
        self.racen = np.zeros(0, dtype=np.float64)
        self.deccen = np.zeros(0, dtype=np.float64)
        self.nfilled = np.zeros(0, dtype=np.float64)
        self.fieldid = np.zeros(0, dtype=np.int32)
        self.nextmjd = np.zeros(0, dtype=np.float64)
        self.cadence = []
        self.observations = []
        self.slots = []
        self.lstObserved = np.zeros(0, dtype=np.int32)
        self.cadencelist = roboscheduler.cadence.CadenceList()
        self._validCadance = None
        self._obsPlan = None
        self._lstPlan = None
        self._lunationPlan = None
        return

    def setPriorities(self):
        scale = [10 if "bhm_rm" in c else 1 for c in self.cadence]
        self.basePriority = self.basePriority * np.array(scale)
        return

    def fromarray(self, fields_array=None):
        self.nfields = len(fields_array)
        self.racen = fields_array['racen']
        self.deccen = fields_array['deccen']
        self.nfilled = fields_array['nfilled']
        self.fieldid = np.arange(self.nfields, dtype=np.int32)
        self.cadence = [c.strip() for c in fields_array['cadence']]
        self.slots = fields_array['slots_exposures']
        self.lstObserved = np.zeros((len(self.slots), 24), dtype=np.int32)
        self.observations = [np.zeros(0, dtype=np.int32)] * self.nfields
        self.icadence = np.zeros(self.nfields, dtype=np.int32)
        self.nextmjd = np.zeros(self.nfields, dtype=np.float64)
        self.basePriority = np.ones(self.nfields) * 200
        self.setPriorities()
        return

    def fromfits(self, filename=None):
        self.fields_fits = fitsio.read(filename)
        self.fromarray(self.fields_fits)
        return

    def add_observations(self, mjd=None, fieldid=None, iobs=None, lst=None):
        self.observations[fieldid] = np.append(self.observations[fieldid],
                                               iobs)
        self.icadence[fieldid] = self.icadence[fieldid] + 1
        cadence = self.cadencelist.cadences[self.cadence[fieldid]]
        if(self.icadence[fieldid] < cadence.nexposures):
            self.nextmjd[fieldid] = (mjd +
                                     cadence.delta_min[self.icadence[fieldid]])
        else:
            self.nextmjd[fieldid] = 100000.
        int_lst = int(np.round(lst/15, 0))
        if int_lst == 24:
            int_lst = 0
        self.lstObserved[fieldid][int_lst] += 1
        return

    @property
    def validCadence(self):
        if self._validCadance is None:
            assert len(self.cadence) > 0, "no field cadences!"
            assert len(self.cadencelist.cadences) > 0, "no cadences in cadencelist!"
            self._validCadance = np.array([c in self.cadencelist.cadences
                                           for c in self.cadence])
        return self._validCadance

    @property
    def obsPlan(self):
        """Slots specify when the field should be observed
        make that information useful.

        Returns:
        -------

        obsPlan : list
            a list containing lst and lunation for each field
        """
        if self._obsPlan is None:
            self._obsPlan = [np.where(s > 0)  for s in self.slots]
        return self._obsPlan

    @property
    def lstPlan(self):
        if self._lstPlan is None:
            # self._lstPlan = np.array([np.mean(p[0]) for p in self.obsPlan])
            self._lstPlan = np.array([np.sum(s, axis=1) for s in self.slots])
        return self._lstPlan

    # @property
    # def lunationPlan(self):
    #     if self._lunationPlan is None:
    #         self._lunationPlan = np.array([np.mean(p[1]) for p in self.obsPlan])
    #     return self._lunationPlan

    def lstWeight(self, lst, fields=None):
        # field id corresponds to indx, so fields is id/indx
        # as is everywhere, but just to keep reminding me...
        assert lst > 0 and lst < 24, "lst must be in hours!"
        diffs = []
        if fields is not None:
            lst_obs = self.lstObserved[fields]
            lst_plan = self.lstPlan[fields]
        else:
            lst_obs = self.lstObserved
            lst_plan = self.lstPlan

        for o, p in zip(lst_obs, lst_plan):
            done = p - o
            elligible = np.where(done > 0)[0]
            if len(elligible) == 0:
                diffs.append(12.)
                continue

            diff = lstDiff(elligible, lst*np.ones(len(elligible)))
            # at the moment we don't care which it lines up with
            diffs.append(np.min(diff))

        return np.array(diffs)

    def toarray(self):
        """Return cadences as a record array

        Returns:
        -------

        fields : ndarray
            information on each field
"""
        maxn = np.array([len(x) for x in self.observations]).max()
        if(maxn == 1):
            maxn = 2
        fields0 = [('fieldid', np.int32),
                   ('racen', np.float64),
                   ('deccen', np.float64),
                   ('cadence', np.dtype('a20')),
                   ('nobservations', np.int32),
                   ('observations', np.int32, maxn)]
        fields = np.zeros(self.nfields, dtype=fields0)
        fields['fieldid'] = self.fieldid
        fields['racen'] = self.racen
        fields['deccen'] = self.deccen
        for indx in np.arange(self.nfields):
            fields['cadence'][indx] = self.cadence[indx]
            fields['nobservations'][indx] = len(self.observations[indx])
            if(fields['nobservations'][indx] > 0):
                fields['observations'][indx][0:fields['nobservations'][indx]] = self.observations[indx]
        return(fields)


def lstDiffSingle(a, b):
    """Intelligently find difference in 2 lsts

    Parameters:
    -----------

    a : np.float32, float
        first lst in hours
    b : np.float32, float
        second lst in hours

    Returns:
    --------

    diff: float
        the absolute difference

    Comments:
    ---------

    """

    if a < b:
        return min(b - a, (a + 24) - b)

    else:
        # a must be bigger
        return min(a - b, (b + 24) - a)


def lstDiff(a, b):
    """wrap lst math to handle arrays

    Parameters:
    -----------

    a : ndarray or list of float32
        first lst in hours
    b : ndarray or list of float32
        second lst in hours

    Returns:
    --------

    diff: ndarray of float32
        the absolute differences

    Comments:
    ---------

    """

    assert len(a) == len(b), "can't compare arrays of different size!"

    return np.array([lstDiffSingle(i, j) for i, j in zip(a, b)])
