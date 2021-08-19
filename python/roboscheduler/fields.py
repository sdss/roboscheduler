import numpy as np
import fitsio
import roboscheduler.cadence
import sys

try:
    import sdssdb.peewee.sdss5db.opsdb as opsdb
    import sdssdb.peewee.sdss5db.targetdb as targetdb
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
    """List of fields. Initialized empty; attributes will be set
    via fromfits() or fromdb().

    Parameters:
    ----------

    plan : string
        robostrategy version, used for db access

    observatory : string
        APO or LCO, used for db access

    Attributes:
    ----------

    nfields : np.int32
        number of fields

    racen : ndarray of np.float64
        right ascension of each design (J2000 deg)

    deccen : ndarray of np.float64
        declination of each design (J2000 deg)

    field_id : ndarray of np.int32
        id for each field

    cadence : list of string
        cadence for each field

    observations : list of ndarray of np.int32
        for each field, indices of observations


    slots : ndarray of np.int32
        Nx2x24, corresponding to robostrategy planned LST bright/dark plan

    lstObserved : ndarray of np.int32
        Nx24 observations in each lst per field;
        respecting robostrategy plan leads to an efficient survey

    cadencelist : robostrategy.cadence.CadenceList singleton
        the cadenceList singleton, for easy access

    """
    def __init__(self, plan=None, observatory=None):
        self.plan = plan
        self.observatory = observatory
        self.nfields = 0
        self.racen = np.zeros(0, dtype=np.float64)
        self.deccen = np.zeros(0, dtype=np.float64)
        self.nfilled = np.zeros(0, dtype=np.float64)
        self.field_id = np.zeros(0, dtype=np.int32)
        # self.nextmjd = np.zeros(0, dtype=np.float64)
        self.epoch_idx = np.zeros(0, dtype=np.int32)
        self.notDone = np.ones(0, dtype=np.bool)  # true is not done to skip the invert elsewhere
        self.cadence = []
        self.observations = []
        self.slots = []
        self.flag = []
        self.lstObserved = np.zeros(0, dtype=np.int32)
        self.cadencelist = roboscheduler.cadence.CadenceList()
        self._validCadance = None
        self._obsPlan = None
        self._lstPlan = None
        self._lunationPlan = None
        self._hist = None
        self._database = _database
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
        # self.nextmjd = np.zeros(self.nfields, dtype=np.float64)
        self.epoch_idx = np.zeros(self.nfields, dtype=np.float64)
        self.basePriority = np.ones(self.nfields) * 200
        self.notDone = np.ones(self.nfields, dtype=np.bool)
        if "flag" in fields_array.dtype.names:
            self.flag = fields_array["flag"]
        else:
            self.flag = np.zeros(self.nfields)
        self.setPriorities()
        return

    def fromfits(self, filename=None):
        """Load a fits file into this Fields object,
           likely an RS-Assignments file.

        Parameters:
        ----------

        filename : string
            the full path to the fits file to be loaded
        """
        self._database = False
        fields_model = [('field_id', np.int32),
                        ('racen', np.float64),
                        ('deccen', np.float64),
                        ('nfilled', np.int32),
                        ('flag', np.int32),
                        ('slots_exposures', np.int32, (24, 2)),
                        ('cadence', np.dtype('a20'))]

        fits_dat = fitsio.read(filename)

        self.fields_fits = np.zeros(len(fits_dat["fieldid"]), dtype=fields_model)

        self.fields_fits["field_id"] = fits_dat["fieldid"]
        self.fields_fits["racen"] = fits_dat["racen"]
        self.fields_fits["deccen"] = fits_dat["deccen"]
        self.fields_fits["nfilled"] = fits_dat["nfilled"]
        self.fields_fits["slots_exposures"] = fits_dat["slots_exposures"]
        self.fields_fits["cadence"] = fits_dat["cadence"]

        self.fromarray(self.fields_fits)
        return

    # def add_observations(self, mjd=None, fieldidx=None, iobs=None,
    #                      lst=None, epoch_idx=None):
    #     # """Add an observation, used in simulations.

    #     # Parameters:
    #     # ----------

    #     # mjd : float
    #     #     MJD of the observation
    #     # fieldidx : integer
    #     #     index into Fields of the field observed
    #     # iobs : integer
    #     #     index into Scheduler.observations
    #     # lst : float
    #     #     lst at time of observation
    #     # epoch_idx : integer
    #     #     which epoch was being observed
    #     # """
    #     # self._hist[self.field_id[fieldidx]].append(mjd)

    #     self.observations[fieldidx] = np.append(self.observations[fieldidx], iobs)
    #     self.icadence[fieldidx] = epoch_idx
    #     cadence = self.cadencelist.cadences[self.cadence[fieldidx]]
    #     if(self.icadence[fieldidx] < cadence.nepochs):
    #         self.nextmjd[fieldidx] = (mjd +
    #                                   cadence.delta_min[self.icadence[fieldidx]])
    #     else:
    #         self.nextmjd[fieldidx] = 100000.
    #         self.icadence[fieldidx] = cadence.nepochs - 1
    #     int_lst = int(np.round(lst/15, 0))
    #     if int_lst == 24:
    #         int_lst = 0
    #     self.lstObserved[fieldidx][int_lst] += 1
    #     return

    @property
    def validCadence(self):
        if self._validCadance is None:
            assert len(self.cadence) > 0, "no field cadences!"
            assert len(self.cadencelist.cadences) > 0, "no cadences in cadencelist!"
            self._validCadance = np.array([c in self.cadencelist.cadences
                                           for c in self.cadence])
        return self._validCadance

    @property
    def lstPlan(self):
        if self._lstPlan is None:
            self._lstPlan = np.array([np.sum(s, axis=1) for s in self.slots])
        return self._lstPlan

    @property
    def hist(self):
        # """The observation history of the field. Unlike Fields.observations,
        # this property will be used on sky.

        # Attributes:
        # ----------

        # _hist : dict of lists
        #     A dictionary with field_id keys containing lists of obs mjds for
        #     each field.
        # """
        if self._hist is None:
            self._hist = {f: list() for f in self.field_id}
            if self._database:
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

            for i in range(len(self.field_id)):
                self.checkCompletion(i)

        return self._hist

    def checkCompletion(self, fieldidx):
        """evaluate field completion
        """
        if not self.notDone[fieldidx]:
            # already done
            return
        field_id = self.field_id[fieldidx]
        cadence = self.cadencelist.cadences[self.cadence[fieldidx]]
        tol = cadence.max_length[int(self.epoch_idx[fieldidx])]
        obs_epochs, begin_last_epoch = epochs_completed(self.hist[field_id],
                                                        tolerance=tol)
        self.epoch_idx[fieldidx] = int(obs_epochs)
        if obs_epochs >= cadence.nepochs:
            self.notDone[fieldidx] = False

    def completeDesign(self, field_idx, mjd, lst, iobs):
        """Maybe just for sims but need a way to add mjds to _hist[field_id]
        """
        self.observations[field_idx] = np.append(self.observations[field_idx], iobs)

        field_id = self.field_id[field_idx]

        self._hist[field_id].append(mjd)

        int_lst = int(np.round(lst/15, 0))
        if int_lst == 24:
            int_lst = 0
        self.lstObserved[field_idx][int_lst] += 1

        self.checkCompletion(field_idx)

    def fromdb(self, version=None):
        """Load this Fields object with fields from the targetdb

        Parameters:
        ----------

        version : string
            db version to grab, if Fields.plan is not set. If passed
            Fields.plan will be reset to version
        """
        if(self._database is False):
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
                        ('flag', np.int32),
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
        flags = list()

        for field in dbfields:
            fieldid.append(field.field_id)
            racen.append(field.racen)
            deccen.append(field.deccen)
            slots_exposures.append(field.slots_exposures)
            cadence.append(field.cadence.label)
            if len(field.priority) > 0:
                if field.priority[0].label == "top":
                    flags.append(1)
                elif field.priority[0].label == "disabled":
                    flags.append(-1)
            else:
                flags.append(0)

        fields = np.zeros(len(dbfields), dtype=fields_model)

        fields["field_id"] = fieldid
        fields["racen"] = racen
        fields["deccen"] = deccen
        fields["flag"] = flags
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
        fields['fieldid'] = self.field_id
        fields['racen'] = self.racen
        fields['deccen'] = self.deccen
        for indx in np.arange(self.nfields):
            fields['cadence'][indx] = self.cadence[indx]
            fields['nobservations'][indx] = len(self.observations[indx])
            if(fields['nobservations'][indx] > 0):
                fields['observations'][indx][0:fields['nobservations'][indx]] = self.observations[indx]
        return(fields)

    def getidx(self, fieldid):
        # return idx into array of fieldid
        return int(np.where(self.field_id == fieldid)[0])


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


def epochs_completed(mjd_past, tolerance=240):
    """Calculate # of observed epochs, allowing
    for more exposures than planned.

    tolerance is in minutes
    """
    if len(mjd_past) == 0:
        return 0, 0

    try:
        len(tolerance)
    except TypeError:
        tolerance = np.array([tolerance])

    epoch_idx = 0

    tolerance = tolerance / 60. / 24.
    begin_last_epoch = mjd_past[0]

    obs_epochs = 1
    prev = begin_last_epoch
    for m in mjd_past:
        delta = m - prev
        if delta < tolerance[epoch_idx]:
            continue
        else:
            obs_epochs += 1
            begin_last_epoch = m
            epoch_idx += 1
            if epoch_idx > len(tolerance) - 1:
                epoch_idx -= 1
        prev = m

    return obs_epochs, begin_last_epoch
