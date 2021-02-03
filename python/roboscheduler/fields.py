import numpy as np
import fitsio
import roboscheduler.cadence
import sys

try:
    import sdssdb.peewee.sdss5db.targetdb as targetdb
    import sdssdb.peewee.sdss5db.opsdb as opsdb
    _database = True
except:
    _database = False


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
    def __init__(self, plan=None, observatory=None):
        self.plan = plan
        self.observatory = observatory
        self.nfields = 0
        self.racen = np.zeros(0, dtype=np.float64)
        self.deccen = np.zeros(0, dtype=np.float64)
        self.nfilled = np.zeros(0, dtype=np.float64)
        self.field_id = np.zeros(0, dtype=np.int32)
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
        self._hist = None
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
        self.field_id = fields_array['field_id']
        self.cadence = [c.strip().decode() for c in fields_array['cadence']]
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

    def add_observations(self, mjd=None, fieldidx=None, iobs=None,
                         lst=None, epoch_idx=None):
        self.observations[fieldidx] = np.append(self.observations[fieldidx], iobs)
        self.icadence[fieldidx] = epoch_idx
        cadence = self.cadencelist.cadences[self.cadence[fieldidx]]
        if(self.icadence[fieldidx] < cadence.nepochs):
            self.nextmjd[fieldidx] = (mjd +
                                      cadence.delta_min[self.icadence[fieldidx]])
        else:
            self.nextmjd[fieldidx] = 100000.
            self.icadence[fieldidx] = cadence.nepochs - 1
        int_lst = int(np.round(lst/15, 0))
        if int_lst == 24:
            int_lst = 0
        self.lstObserved[fieldidx][int_lst] += 1
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
            self._obsPlan = [np.where(s > 0) for s in self.slots]
        return self._obsPlan

    @property
    def lstPlan(self):
        if self._lstPlan is None:
            self._lstPlan = np.array([np.sum(s, axis=1) for s in self.slots])
        return self._lstPlan

    # @property
    # def lunationPlan(self):
    #     if self._lunationPlan is None:
    #         self._lunationPlan = np.array([np.mean(p[1]) for p in self.obsPlan])
    #     return self._lunationPlan

    @property
    def hist(self):
        if self._hist is None:
            self._hist = {f: list() for f in self.field_id}
            if _database:
                versionDB = targetdb.Version()
                ver = versionDB.get(plan=self.plan)

                obsDB = targetdb.Observatory()
                obs = obsDB.get(label=self.observatory.upper())

                Field = targetdb.Field
                Design = targetdb.Design
                Status = opsdb.CompletionStatus
                d2s = opsdb.DesignToStatus
                done = Status.get(label="done")

                dbfields = Field.select(d2s.mjd, Field.field_id)\
                                .join(Design)\
                                .join(d2s, on=(Design.pk == d2s.design_pk))\
                                .where((Field.version == ver) &
                                       (Field.observatory == obs),
                                       (d2s.status == done)).dicts()

                for d in dbfields:
                    self._hist[d["field_id"]].append(d["mjd"])

        return self._hist

    def fromdb(self, version=None):
        """Extract cadences into the targetdb
        """
        if(_database is False):
            print("No database available.")
            return()

        if version is None:
            version = self.plan
        else:
            self.plan = version
        assert version is not None, "must specify version!"

        fields_model = [('field_id', np.int32),
                        ('racen', np.float64),
                        ('deccen', np.float64),
                        ('nfilled', np.int32),
                        ('slots_exposures', np.int32, (24, 2)),
                        ('cadence', np.dtype('a20'))]

        versionDB = targetdb.Version()
        ver = versionDB.get(plan=version)

        obsDB = targetdb.Observatory()
        obs = obsDB.get(label=self.observatory.upper())

        Field = targetdb.Field
        dbfields = Field.select().where(Field.version == ver,
                                        Field.observatory == obs)

        fieldid = list()
        racen = list()
        deccen = list()
        slots_exposures = list()
        cadence = list()

        for field in dbfields:
            fieldid.append(field.field_id)
            racen.append(field.racen)
            deccen.append(field.deccen)
            slots_exposures.append(field.slots_exposures)
            cadence.append(field.cadence.label)

        fields = np.zeros(len(dbfields), dtype=fields_model)

        fields["field_id"] = fieldid
        fields["racen"] = racen
        fields["deccen"] = deccen
        fields["slots_exposures"] = slots_exposures
        fields["cadence"] = cadence

        self.fromarray(fields_array=fields)

        self.cadencelist.fromdb()

    def lstWeight(self, lst, field_idx=None):
        # field id corresponds to indx, so fields is id/indx
        # as is everywhere, but just to keep reminding me...
        assert lst > 0 and lst < 24, "lst must be in hours!"
        diffs = []
        if field_idx is not None:
            lst_obs = self.lstObserved[field_idx]
            lst_plan = self.lstPlan[field_idx]
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
