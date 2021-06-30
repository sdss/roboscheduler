import numpy as np

"""Observations module class.

Dependencies:

 numpy
 scipy

"""


def wrapHA(ha):
    """Force ha into range -180 to 180. Just to be sure.
    """
    if -180 < ha <= 180:
        return ha
    while ha > 180:
        ha = ha - 360.
    while ha <= -180:
        ha = ha + 360.
    assert -180 < ha <= 180, "ha = {:.2f}".format(ha)
    return ha


class Observations(object):
    """Observations class.

    Parameters:
    ----------

    observatory : string
       Observatory for fields ('apo' or 'lco'; default 'apo')

    Attributes:
    ----------

    nobservations : np.int32
        number of observations

    fieldid : ndarray of np.int32
        id of each field for observations

    mjd : ndarray of np.float64
        MJD of observation (days)

    duration : ndarray of np.float64
        duration of observation (days)

    sn2 : ndarray of np.float32
        duration of observation (days)

    airmass : ndarray of np.float32
        airmass of observation

    skybrightness : ndarray of np.float32
        sky brightness of observation (Moon illumination fraction, or zero if
        below horizon, or one if in bright twilight)

    lst : ndarray of np.float32
        LST of observation

    Methods:
    -------

    add() : add an observation of a field
    toarray() : Return ndarray of field properties

"""
    def __init__(self, observatory='apo'):
        self.nobservations = np.int32(0)
        self.observatory = observatory
        self.fieldid = np.zeros(0, dtype=np.int32)
        self.racen = np.zeros(0, dtype=np.float64)
        self.deccen = np.zeros(0, dtype=np.float64)
        self.cadence = np.zeros(0, dtype=np.dtype('a20'))
        self.nfilled = np.zeros(0, dtype=np.int32)
        self.mjd = np.zeros(0, dtype=np.float64)
        self.duration = np.zeros(0, dtype=np.float64)
        self.apgSN2 = np.zeros(0, dtype=np.float32)
        self.rSN2 = np.zeros(0, dtype=np.float32)
        self.bSN2 = np.zeros(0, dtype=np.float32)
        self.airmass = np.zeros(0, dtype=np.float32)
        self.skybrightness = np.zeros(0, dtype=np.float32)
        self.lst = np.zeros(0, dtype=np.float32)
        self.ha = np.zeros(0, dtype=np.float64)
        self.nexp_cumul = np.zeros(0, dtype=np.int32)
        return

    def add(self, fieldid=None, mjd=None, duration=None, apgSN2=None,
            rSN2=None, bSN2=None, skybrightness=None, airmass=None,
            lst=None, racen=None, deccen=None, cadence=None, nfilled=None,
            nexp_cumul=None):
        self.fieldid = np.append(self.fieldid,
                                 np.array([np.float64(fieldid)]))
        self.mjd = np.append(self.mjd,
                             np.array([np.float64(mjd)]))
        self.duration = np.append(self.duration,
                                  np.array([np.float64(duration)]))
        self.apgSN2 = np.append(self.apgSN2,
                                np.array([np.float32(apgSN2)]))
        self.rSN2 = np.append(self.rSN2,
                              np.array([np.float32(rSN2)]))
        self.bSN2 = np.append(self.bSN2,
                              np.array([np.float32(bSN2)]))
        self.skybrightness = np.append(self.skybrightness,
                                       np.array([np.float32(skybrightness)]))
        self.airmass = np.append(self.airmass,
                                 np.array([np.float32(airmass)]))
        self.lst = np.append(self.lst,
                             np.array([np.float32(lst)]))
        self.nobservations = len(self.fieldid)

        self.racen = np.append(self.racen,
                               np.array([np.float64(racen)]))
        self.deccen = np.append(self.deccen,
                                np.array([np.float64(deccen)]))
        self.cadence = np.append(self.cadence,
                                 np.array([cadence]))
        self.nfilled = np.append(self.nfilled,
                                 np.array([np.int32(nfilled)]))

        ha = wrapHA(lst - racen)
        self.ha = np.append(self.ha,
                            np.array([np.float64(ha)]))
        self.nexp_cumul = np.append(self.nexp_cumul,
                                    np.array([np.int32(nexp_cumul)]))

        return(self.nobservations - 1)

    def forfield(self, mjd=None, fieldid=None):
        indx = np.where((self.mjd <= mjd) &
                        (self.fieldid == fieldid))[0]
        return(self.toarray(indx=indx))

    def toarray(self, indx=None):
        """Return observations as a record array

        Parameters:
        ----------

        indx : ndarray of np.int32
            indices of observations to return (default to all)

        Returns:
        -------

        observations : record array
            observation information
        """
        obs0 = [('fieldid', np.int32),
                ('mjd', np.float64),
                ('duration', np.float64),
                ('apgSN2', np.float32),
                ('rSN2', np.float32),
                ('bSN2', np.float32),
                ('airmass', np.float32),
                ('skybrightness', np.float32),
                ('lst', np.float32),
                ('racen', np.float64),
                ('deccen', np.float64),
                ('cadence', np.dtype('a20')),
                ('nfilled', np.int32),
                ('ha', np.float32),
                ('nexp_cumul', np.int32)]
        if(indx is None):
            indx = np.arange(self.nobservations)
        nobs = len(indx)
        obs = np.zeros(nobs, dtype=obs0)
        if(nobs > 0):
            obs['fieldid'] = self.fieldid[indx]
            obs['mjd'] = self.mjd[indx]
            obs['duration'] = self.duration[indx]
            obs['apgSN2'] = self.apgSN2[indx]
            obs['rSN2'] = self.rSN2[indx]
            obs['bSN2'] = self.bSN2[indx]
            obs['airmass'] = self.airmass[indx]
            obs['skybrightness'] = self.skybrightness[indx]
            obs['lst'] = self.lst[indx]
            obs['racen'] = self.racen[indx]
            obs['deccen'] = self.deccen[indx]
            obs['cadence'] = self.cadence[indx]
            obs['nfilled'] = self.nfilled[indx]
            obs['ha'] = self.ha[indx]
            obs['nexp_cumul'] = self.nexp_cumul[indx]
        return(obs)
