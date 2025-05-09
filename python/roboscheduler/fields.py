import os

import numpy as np
import fitsio
import roboscheduler.cadence

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

    scheduler: roboscheduler.scheduler.Scheduler or None
        Scheduler object for on sky math

    Attributes:
    ----------

    nfields : np.int32
        number of fields

    racen : ndarray of np.float64
        right ascension of each design (J2000 deg)

    deccen : ndarray of np.float64
        declination of each design (J2000 deg)

    pk : ndarray of np.int32
        db primary key for each field, changed from field_id to allow
        multiple fields with the same id (and therefore position on sky)
        to be scheduled together

    cadence : list of string
        cadence for each field

    observations : list of ndarray of np.int32
        for each field, indices of observations


    slots : ndarray of np.int32
        Nx2x24, corresponding to robostrategy planned LST bright/dark plan

    lstObserved : ndarray of np.int32
        Nx24x2 observations in each lst per field;
        respecting robostrategy plan leads to an efficient survey

    cadencelist : robostrategy.cadence.CadenceList singleton
        the cadenceList singleton, for easy access

    """
    def __init__(self, plan=None, observatory=None, scheduler=None):
        self.plan = plan
        self.observatory = observatory
        self.nfields = 0
        self.racen = np.zeros(0, dtype=np.float64)
        self.deccen = np.zeros(0, dtype=np.float64)
        self.nfilled = np.zeros(0, dtype=np.float64)
        self.pk = np.zeros(0, dtype=np.int32)
        # self.nextmjd = np.zeros(0, dtype=np.float64)
        self.epoch_idx = np.zeros(0, dtype=np.int32)
        self.notDone = np.ones(0, dtype=np.bool_)  # true is not done to skip the invert elsewhere
        self.cadence = []
        self.observations = []
        self.slots = []
        self.flag = []
        self.lstObserved = np.zeros(0, dtype=np.int32)
        self.cadencelist = roboscheduler.cadence.CadenceList(observatory=observatory)
        self.scheduler = scheduler
        self._validCadance = None
        self._obsPlan = None
        self._lstPlan = None
        self._lunationPlan = None
        self._hist = None
        self._database = _database
        self._designs = None
        return

    def fromarray(self, fields_array=None, designList=None):
        self.nfields = len(fields_array)
        self.racen = fields_array['racen']
        self.deccen = fields_array['deccen']
        self.nfilled = fields_array['nfilled']
        self.field_id = fields_array['field_id']
        self.pk = fields_array['pk']
        self.cadence = [c.strip().decode() for c in fields_array['cadence']]
        self.slots = fields_array['slots_exposures']
        self.lstObserved = fields_array["lstObserved"]
        self.observations = [np.zeros(0, dtype=np.int32)] * self.nfields
        self.icadence = np.zeros(self.nfields, dtype=np.int32)
        # self.nextmjd = np.zeros(self.nfields, dtype=np.float64)
        self.epoch_idx = np.zeros(self.nfields, dtype=np.float64)
        self.original_exposures_done = fields_array["original_exposures_done"]
        if "base_priority" in fields_array.dtype.names:
            self.basePriority = fields_array["base_priority"]
        else:
            self.basePriority = np.ones(self.nfields)
        self.notDone = np.ones(self.nfields, dtype=np.bool_)
        if "flag" in fields_array.dtype.names:
            self.flag = fields_array["flag"]
        else:
            self.flag = np.zeros(self.nfields)
        # self.setPriorities()
        if designList:
            assert len(designList) == len(self.pk), "designList must match fields"
            self._designs = designList
        return

    def createDummyDesigns(self):
        """Create a list of arrays of dummy design ids
        """
        runningNumber = 1e6
        self._designs = list()
        for n in self.nfilled:
            self._designs.append(np.arange(n) + runningNumber)
            runningNumber += n
            assert runningNumber not in self._designs[-1], "bad math"

    @property
    def designs(self):
        if self._designs is None:
            self.createDummyDesigns()
        return self._designs

    def fromfits(self, filename=None, priority_fields={}):
        """Load a fits file into this Fields object,
           likely an RS-Allocation file.

        Parameters:
        ----------

        filename : string
            the full path to the fits file to be loaded
        """
        fits_dat = fitsio.read(filename)

        try:
            len_exposures = fits_dat["original_exposures_done"].shape[1]
        except ValueError:
            len_exposures = 1

        self._database = False
        fields_model = [('pk', np.int32),
                        ('field_id', np.int32),
                        ('racen', np.float64),
                        ('deccen', np.float64),
                        ('nfilled', np.int32),
                        ('flag', np.int32),
                        ('slots_exposures', np.int32, (24, 2)),
                        ('lstObserved', np.int32, (24, 2)),
                        ('original_exposures_done', np.int32, (len_exposures)),
                        ('cadence', np.dtype('a40')),
                        ('base_priority', np.int32)]

        self.fields_fits = np.zeros(len(fits_dat["fieldid"]), dtype=fields_model)

        self.fields_fits["pk"] = np.arange(len(fits_dat["fieldid"]))
        self.fields_fits["field_id"] = fits_dat["fieldid"]
        self.fields_fits["racen"] = fits_dat["racen"]
        self.fields_fits["deccen"] = fits_dat["deccen"]
        if "nfilled" in fits_dat.dtype.names:
            self.fields_fits["nfilled"] = fits_dat["nfilled"]
        elif "nallocated_full" in fits_dat.dtype.names:
            self.fields_fits["nfilled"] = fits_dat["nallocated_full"]
        else:
            print("WARN: strange rsAllocation format, estimating nfilled")
            self.fields_fits["nfilled"] =\
                  np.round(np.sum(np.sum(fits_dat["slots_exposures"], axis=1),
                         axis=1)).astype(int)
        self.fields_fits["slots_exposures"] = fits_dat["slots_exposures"]
        self.fields_fits["cadence"] = fits_dat["cadence"]
        if len_exposures > 1:
            self.fields_fits["original_exposures_done"] = fits_dat["original_exposures_done"]

        self.fields_fits["base_priority"] = np.ones(len(fits_dat["fieldid"]))

        for i, f in enumerate(self.fields_fits):
            if f["field_id"] in priority_fields:
                if str(f["cadence"].decode()) == str(priority_fields[f["field_id"]]["cadence"]):
                    self.fields_fits["base_priority"][i] += priority_fields[f["field_id"]]["priority"]

        w_cvz = np.where(["100x8" in c for c in fits_dat["cadence"]])

        self.fields_fits["flag"][w_cvz] = -1

        self.fromarray(self.fields_fits)
        self.createDummyDesigns()
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
        #     A dictionary with pk keys containing lists of obs mjds for
        #     each field.
        # """
        if self._hist is None:
            self._hist = {f: list() for f in self.pk}
            if self._database:
                versionDB = targetdb.Version()
                ver = versionDB.get(plan=self.plan)

                obsDB = targetdb.Observatory()
                obs = obsDB.get(label=self.observatory.upper())

                Field = targetdb.Field
                Design = targetdb.Design
                d2f = targetdb.DesignToField
                Status = opsdb.CompletionStatus
                d2s = opsdb.DesignToStatus
                done = Status.get(label="done")

                dbfields = Field.select(d2s.mjd, Field.pk)\
                                .join(d2f, on=(Field.pk == d2f.field_pk))\
                                .join(Design, on=(Design.design_id == d2f.design_id))\
                                .join(d2s, on=(Design.design_id == d2s.design_id))\
                                .where((Field.version == ver) &
                                       (Field.observatory == obs),
                                       (d2s.status == done)).dicts()

                for d in dbfields:
                    self._hist[d["pk"]].append(d["mjd"])

        return self._hist

    def checkCompletion(self, fieldidx):
        """evaluate field completion
        """
        if (not self.notDone[fieldidx]
            or not self.validCadence[fieldidx]):
            # already done
            return
        pk = self.pk[fieldidx]
        cadence = self.cadencelist.cadences[self.cadence[fieldidx]]
        tol = cadence.max_length[int(self.epoch_idx[fieldidx])]
        obs_epochs, begin_last_epoch = epochs_completed(self.hist[pk],
                                                        tolerance=tol,
                                                        nexp=cadence.nexp[0])
        self.epoch_idx[fieldidx] = int(obs_epochs)
        if obs_epochs >= cadence.nepochs:
            self.notDone[fieldidx] = False

    def completeDesign(self, field_idx, mjd, lst, iobs, dark=False):
        """Maybe just for sims but need a way to add mjds to _hist[pk]
        """
        self.observations[field_idx] = np.append(self.observations[field_idx], iobs)

        pk = self.pk[field_idx]

        self._hist[pk].append(mjd)

        if dark:
            row = 0
        else:
            row = 1

        int_lst = int(np.round(lst/15, 0))
        if int_lst == 24:
            int_lst = 0
        self.lstObserved[field_idx][int_lst, row] += 1

        self.checkCompletion(field_idx)

    def fromdb(self, version=None, priorities={}):
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

        fields_model = [('pk', np.int32),
                        ('field_id', np.int32),
                        ('racen', np.float64),
                        ('deccen', np.float64),
                        ('nfilled', np.int32),
                        ('flag', np.int32),
                        ('base_priority', np.int32),
                        ('slots_exposures', np.int32, (24, 2)),
                        ('lstObserved', np.int32, (24, 2)),
                        ('original_exposures_done', np.int32, (1)),
                        ('cadence', np.dtype('a40'))]

        versionDB = targetdb.Version()
        ver = versionDB.get(plan=version)

        obsDB = targetdb.Observatory()
        obs = obsDB.get(label=self.observatory.upper())

        Field = targetdb.Field
        dbfields = Field.select().where(Field.version == ver,
                                        Field.observatory == obs)

        field_pri_ver = os.getenv('FIELD_PRI_VER')
        if field_pri_ver is None:
            field_pri_ver = ver

        pri_ver = opsdb.PriorityVersion.get(label=field_pri_ver)

        bp = opsdb.BasePriority
        prioritizedFields = bp.select().where(bp.version==pri_ver).dicts()
        priorityDict = {p["field"]: p["priority"] for p in prioritizedFields}

        lstHist = opsdb.LstHist

        lstDone = lstHist.select(lstHist.field, lstHist.lst_counts)\
                         .join(Field)\
                         .where(Field.version == ver).dicts()

        lstDict = {l["field"]: l["lst_counts"] for l in lstDone}

        pk = list()
        field_id = list()
        racen = list()
        deccen = list()
        slots_exposures = list()
        lstObserved = list()
        cadence = list()
        flags = list()
        base_priority = list()

        for field in dbfields:
            pk.append(field.pk)
            field_id.append(field.field_id)
            racen.append(field.racen)
            deccen.append(field.deccen)
            slots_exposures.append(field.slots_exposures)
            lstObserved.append(lstDict[field.pk])
            cadence.append(field.cadence.label)
            if len(field.priority) > 0:
                if field.priority[0].label == "top":
                    flags.append(1)
                elif field.priority[0].label == "disabled":
                    flags.append(-1)
            else:
                flags.append(0)
            if field.pk in priorityDict:
                base_priority.append(priorityDict[field.pk])
            else:
                base_priority.append(1)

        fields = np.zeros(len(dbfields), dtype=fields_model)

        fields["pk"] = pk
        fields["field_id"] = field_id
        fields["racen"] = racen
        fields["deccen"] = deccen
        fields["flag"] = flags
        fields["base_priority"] = base_priority
        fields["slots_exposures"] = slots_exposures
        fields["lstObserved"] = lstObserved
        fields["cadence"] = cadence

        # we're resetting field hist, so need to re-cache
        self._hist = None

        self.fromarray(fields_array=fields)

        self.cadencelist.fromdb(use_label_root=False, version="v2",
                                priorities=priorities)

    def lstWeight(self, lst, field_idx=None, dark=False):
        # field id corresponds to indx, so fields is id/indx
        # as is everywhere, but just to keep reminding me...
        assert lst > 0 and lst < 24, "lst must be in hours!"
        if dark:
            row = 0
        else:
            row = 1
        if field_idx is not None:
            # lst_obs = self.lstObserved[field_idx]
            # lst_plan = self.lstPlan[field_idx]
            lst_plan = self.slots[field_idx][:, :, row]
            lst_obs = self.lstObserved[field_idx][:, :, row]
        # I don't think this is used, let it break if so
        else:
            # assert False, "you called the wrong LST weight!"
            lst_plan = self.slots[:, :, row]
            lst_obs = self.lstObserved[:, :, row]

        diffs = []
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
        fields0 = [('pk', np.int32),
                   ('fieldid', np.int32),
                   ('racen', np.float64),
                   ('deccen', np.float64),
                   ('cadence', np.dtype('a40')),
                   ('nobservations', np.int32),
                   ('nfilled', np.int32),
                   ('observations', np.int32, maxn),
                   ('base_priority', np.int32)]
        fields = np.zeros(self.nfields, dtype=fields0)
        fields['pk'] = self.pk
        fields['fieldid'] = self.field_id
        fields['racen'] = self.racen
        fields['nfilled'] = self.nfilled
        fields['deccen'] = self.deccen
        fields['base_priority'] = self.basePriority
        for indx in np.arange(self.nfields):
            fields['cadence'][indx] = self.cadence[indx]
            fields['nobservations'][indx] = len(self.observations[indx])
            if(fields['nobservations'][indx] > 0):
                fields['observations'][indx][0:fields['nobservations'][indx]] = self.observations[indx]
        return(fields)

    def getidx(self, fieldid):
        # return idx into array of fieldid
        return int(np.where(self.pk == fieldid)[0])


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


def epochs_completed(mjd_past, tolerance=0.5, nexp=1):
    """Calculate # of observed epochs, allowing
    for more exposures than planned.

    tolerance HAS CHANGED TO DAYS
    """
    if len(mjd_past) < nexp:
        return 0, 0

    epoch_idx = 0

    tolerance = tolerance
    begin_last_epoch = mjd_past[0]

    obs_epochs = 1
    prev = begin_last_epoch
    for m in mjd_past:
        delta = m - prev
        if delta < tolerance:
            continue
        else:
            obs_epochs += 1
            begin_last_epoch = m
            epoch_idx += 1
        prev = m

    return obs_epochs, begin_last_epoch
