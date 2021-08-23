import os
import numpy as np
import scipy.optimize as optimize
import PyAstronomy.pyasl as pyasl
import astropy.units as units
import astropy.time as atime
import pydl.pydlutils.yanny as yanny
import roboscheduler.fields
import roboscheduler.observations
import roboscheduler.cadence
from roboscheduler.moonphase import moonphase2
from roboscheduler.sunpos2 import sunpos2
from roboscheduler.ks91 import KS91_deltaV
from roboscheduler.fields import epochs_completed


"""Scheduler module class.

Dependencies:

 numpy
 scipy
 PyAstronomy
 astropy
 pydl

"""


def dateandtime2mjd(date=None, time='12:00', to_tai=7):
    """Utility to calculate an MJD"""
    if((type(date) is list) | (type(date) is np.ndarray)):
        isotimes = ["{date} {time}".format(date=cdate, time=ctime)
                    for cdate, ctime in zip(date, time)]
    else:
        isotimes = "{date} {time}".format(date=date, time=time)
    times = atime.Time(isotimes, format='iso', scale='tai')
    times = times + np.int32(to_tai) * units.hour
    return(times.mjd)


class SchedulerBase(object):
    """Scheduler base class with generic utilities.

    Parameters:
    ----------

    Attributes:
    ----------

    Methods:
    -------

    ralst2ha(ra=, lst=) : convert RA and LST to hour angle
    hadec2altaz(ha=, dec=, lat=) : convert HA, Dec, latitude to alt, az
    alt2airmass(alt=) : convert altitude to airmass

    """
    def __init__(self):
        return

    def _arrayify(self, quantity=None):
        """Cast quantity as ndarray of numpy.float64"""
        try:
            length = len(quantity)
        except TypeError:
            length = 1
        return np.zeros(length, dtype=np.float64) + quantity

    def _mjd2jd(self, mjd=None):
        """Convert MJD to JD"""
        return (self._arrayify(mjd) + np.float64(2400000.5))

    def ralst2ha(self, ra=None, lst=None):
        """Return HA (degrees) given RA and LST

        Parameters:
        ----------

        ra : np.float64
            right ascension (deg)

        lst : np.float64
            local sidereal time (deg)

        Returns:
        -------

        ha : np.float64
            hour angle (deg)

        """
        ha = (((self._arrayify(lst) - self._arrayify(ra) + 360. + 180.)
               % 360.) - 180.)
        return(ha)

    def hadec2altaz(self, ha=None, dec=None, lat=None):
        """Return (alt, az) (degrees) of given HA and Dec and latitude

        Parameters:
        ----------

        ha : np.float64
            hour angle (deg)

        dec : np.float64
            declination (deg)

        lat : np.float64
            latitude (deg)

        Returns:
        -------

        alt : np.float64
            altitude (deg)

        az : np.float64
            azimuth (deg E of N)
        """
        aha = self._arrayify(ha)
        adec = self._arrayify(dec)
        (alt, az) = pyasl.hadec2altaz(aha, adec,
                                      np.float64(lat) + np.zeros(len(aha)))
        return (alt, az)

    def alt2airmass(self, alt):
        """Return airmass given altitude

        Parameters:
        ----------

        alt : np.float64
            altitude (deg)

        Returns:
        -------

        airmass : np.float64
            airmass (1/sin(altitude))

        """
        airmass = 1. / np.sin(np.pi / 180. * self._arrayify(alt))
        return(airmass)


class Observer(SchedulerBase):
    """Observer class to define different observatories and compute useful 
    quantities (lst, moon location, etc) at each observatory.

    Parameters:
    ----------

    observatory : str
        Name of observatory to use (must be in observatory file)
        (default 'apo')

    observatoryfile : str
        Name of Yanny-format observatory file to read
        (default $ROBOSCHEDULER_DIR/data/observatories.par)

    Attributes:
    ----------

    observatory : str
        Name of observatory

    latitude : numpy.float64
        Latitude of observatory

    longitude : numpy.float64
        Longitude (E of Greenwich) of observatory

    """
    def __init__(self, observatory='apo', observatoryfile=None,
                 dark_twilight=-15., bright_twilight=-8.):
        """Create Observer object"""
        super().__init__()
        self.observatory = observatory
        if(observatoryfile is None):
            prod_dir = os.path.abspath(__file__).split("/scheduler.py")[0]
            observatoryfile = os.path.join(prod_dir,
                                           'etc', 'observatories.par')
        self._file = observatoryfile
        self._data = yanny.yanny(self._file)
        observatories = np.array([obs.decode()
                                  for obs in
                                  self._data['OBSERVATORY']['observatory']])
        indx = np.where(observatories == self.observatory)[0]
        self.latitude = self._data['OBSERVATORY']['latitude'][indx]
        self.longitude = self._data['OBSERVATORY']['longitude'][indx]
        self.dark_twilight = np.float32(dark_twilight)
        self.bright_twilight = np.float32(bright_twilight)
        return

    def lst(self, mjd=None):
        """Return LST (degrees) given MJD for observer

        Parameters:
        ----------

        mjd : np.float64
            Modified Julian Day (days)

        Returns:
        -------

        lst : np.float64
            local sidereal time (deg)
        """
        mjds = self._arrayify(mjd)
        lst = (np.float64(15.) *
               pyasl.ct2lst(self._mjd2jd(mjds),
                            np.zeros(len(mjds)) + self.longitude))
        return (lst)

    def sun_radec(self, mjd=None):
        """Return (ra, dec) in deg J2000 for Sun at MJD (days)

        Parameters:
        ----------

        mjd : np.float64
            Modified Julian Day (days)

        Returns:
        -------

        ra : np.float64
            right ascension, J2000 (deg)

        dec : np.float64
            declination, J2000 (deg)
        """
        jd = self._mjd2jd(mjd=self._arrayify(mjd))
        (tmp_jd, ra, dec) = sunpos2(jd)
        return (ra, dec)

    def moon_radec(self, mjd=None):
        """Return (ra, dec) in deg J2000 for Moon at MJD (days)

        Parameters:
        ----------

        mjd : np.float64
            Modified Julian Day (days)

        Returns:
        -------

        ra : np.float64
            right ascension, J2000 (deg)

        dec : np.float64
            declination, J2000 (deg)
        """
        jd = self._mjd2jd(mjd=self._arrayify(mjd))
        ra, dec, dist, geolon, geolat = pyasl.moonpos(jd)
        return (ra, dec)

    def moon_dist(self, mjd=None, ra=None, dec=None, threshold=None):
        """Return distance to the moon in deg at MJD (days)

        Parameters:
        ----------

        mjd : np.float64
            Modified Julian Day (days)

        ra : np.float64
            right ascension, J2000 (deg)

        dec : np.float64
            declination, J2000 (deg)

        threshold : np.float64
            minimum allowed distance (deg); if provided a boolean array
            of the distance check is returned, other the distance array is returned

        Returns:
        -------

        distance : np.float64
            angular distance to the moon (deg)

        far_enough: np.bool
            array of distance checks for input targets,
            returned if threshold is passed
        """
        ra = self._arrayify(ra)
        dec = self._arrayify(dec)

        moonra, moondec = self.moon_radec(mjd)

        moon_dist = np.power((ra - moonra)*np.cos(dec*np.pi / 180), 2)\
                  + np.power((dec - moondec), 2)

        if threshold:
            return moon_dist > np.power(threshold, 2)
        else:
            return np.power(moon_dist, 0.5)

    def radec2altaz(self, mjd=None, ra=None, dec=None):
        """Return (alt, az) for (ra, dec) in deg J2000 at MJD (days)

        Parameters:
        ----------

        mjd : np.float64
            Modified Julian Day (days)

        ra : np.float64
            right ascension, J2000 (deg)

        dec : np.float64
            declination, J2000 (deg)

        Returns:
        -------

        alt : np.float64
            altitude (deg)

        az : np.float64
            azimuth (deg E of N)
        """
        lst = self.lst(mjd=mjd)
        ha = self.ralst2ha(ra=ra, lst=lst)
        (alt, az) = self.hadec2altaz(ha=ha, dec=dec, lat=self.latitude)
        return (alt, az)

    def moon_illumination(self, mjd=None):
        """Return Moon illumination at MJD (days)

        Parameters:
        ----------

        mjd : np.float64
            Modified Julian Day (days)

        Returns:
        -------

        illumination : np.float64
            fraction of Moon illuminated
        """
        jd = self._mjd2jd(mjd=self._arrayify(mjd))
        return (moonphase2(jd))

    def sun_altaz(self, mjd=None):
        """Return (alt, az) for Sun at MJD (days)

        Parameters:
        ----------

        mjd : np.float64
            Modified Julian Day (days)

        Returns:
        -------

        alt : np.float64
            altitude (deg)

        az : np.float64
            azimuth (deg E of N)
        """
        (ra, dec) = self.sun_radec(mjd=mjd)
        (alt, az) = self.radec2altaz(mjd=mjd, ra=ra, dec=dec)
        return (alt, az)

    def moon_altaz(self, mjd=None):
        """Return (alt, az) for Moon at MJD (days)

        Parameters:
        ----------

        mjd : np.float64
            Modified Julian Day (days)

        Returns:
        -------

        alt : np.float64
            altitude (deg)

        az : np.float64
            azimuth (deg E of N)
        """
        (ra, dec) = self.moon_radec(mjd=mjd)
        (alt, az) = self.radec2altaz(mjd=mjd, ra=ra, dec=dec)
        return (alt, az)

    def lunation(self, mjd=None):
        """Return Moon illumination, or zero if Moon at alt<0"""
        (moon_alt, moon_az) = self.moon_altaz(mjd=mjd)
        if(moon_alt < 0):
            return(0.)
        else:
            return(self.moon_illumination(mjd=mjd))

    def skybrightness(self, mjd=None):
        """Return a sky brightness related number"

        Parameters:
        ----------

        mjd : np.float64
            Modified Julian Day (days)

        Returns:
        -------

        skybrightness : np.float32
            sky brightness related number between 0 and 1

        Notes:
        -----

        If the Sun is above Scheduler.dark_twilight, then the
        skybright is one. Otherwise the skybrightness is equal to the
        lunation, which if the Moon is above the horizon, is its
        fractional illumination, and if the Moon is below the horizon,
        is zero.
        """
        (moon_alt, moon_az) = self.moon_altaz(mjd=mjd)
        (sun_alt, sun_az) = self.sun_altaz(mjd=mjd)
        if(sun_alt > self.dark_twilight):
            return(1.)
        else:
            return(self.lunation(mjd=mjd))

    def _twilight_function(self, mjd=None, twilight=-8.):
        """Utility function for root-finding to get twilight times"""
        (alt, az) = self.sun_altaz(mjd=mjd)
        return (alt - twilight)

    def evening_twilight(self, mjd=None, twilight=None):
        """Return MJD (days) of evening twilight for MJD

        Parameters:
        ----------

        mjd : np.int32, int
            Modified Julian Day (days)

        Returns:
        -------

        evening_twilight : np.float64
            time of twilight in MJD (days)
        """
        if twilight is None:
            twilight = self.bright_twilight
        if(np.floor(np.float64(mjd)) != np.float64(mjd)):
            raise ValueError("MJD should be an integer")
        noon_ish = (np.float64(mjd) -
                    self.longitude / 15. / 24. - 0.5)
        midnight_ish = noon_ish + 0.5
        twi = optimize.brenth(self._twilight_function,
                              noon_ish, midnight_ish,
                              args=twilight)
        return(np.float64(twi))

    def morning_twilight(self, mjd=None, twilight=None):
        """Return MJD (days) of morning twilight for MJD

        Parameters:
        ----------

        mjd : np.int32, int
            Modified Julian Day (days)

        Returns:
        -------

        morning_twilight : np.float64
            time of twilight in MJD (days)
        """
        if twilight is None:
            twilight = self.bright_twilight
        if(np.floor(np.float64(mjd)) != np.float64(mjd)):
            raise ValueError("MJD should be an integer")
        midnight_ish = (np.float64(mjd) -
                        self.longitude / 15. / 24.)
        nextnoon_ish = midnight_ish + 0.5
        twi = optimize.brenth(self._twilight_function,
                              midnight_ish, nextnoon_ish,
                              args=twilight)
        return(np.float64(twi))

    def _moon_rise_set(self, mjd=None):
        """Utility function for root-finding to get moon rise/set times"""
        (alt, az) = self.moon_altaz(mjd=mjd)
        return alt

    def moon_rise_set(self, mjd=None):
        """Return MJD (days) of next moon traverse of the horizon
           and whether it was rise or set.

        Parameters:
        ----------

        mjd : np.int32, int
            Modified Julian Day (days)

        Returns:
        -------

        event : np.float64
            time of moon_rise_set in MJD (days)
        """
        alt, az = self.moon_altaz(mjd=mjd)

        if np.isclose(alt, 0, atol=0.1):
            next_brightness = np.float64(self.moon_illumination(mjd=mjd+0.05))
            return mjd, next_brightness

        # by default it should traverse once in the next 12 hrs
        rise = False
        if alt < 0:
            rise = True
            if az < 180:
                # print(f"{float(alt):.1f}, {float(az):.1f} it's gonna rise soon, be safe, don't get the set too")
                guess = 0.4
            else:
                # print(f"{float(alt):.1f}, {float(az):.1f} it set recently")
                guess = 0.7

        elif az > 180:
            # print(f"{float(alt):.1f}, {float(az):.1f} it's on the setting side")
            guess = 0.4
        else:
            # print(f"{float(alt):.1f}, {float(az):.1f} it's on the rising side, be up awhile")
            guess = 0.7

        event = optimize.brenth(self._moon_rise_set,
                                mjd, mjd+guess)

        # try:
        #     event = optimize.brenth(self._moon_rise_set,
        #                         mjd, mjd+guess)
        # except:
        #     guesses = np.arange(0, 1, 0.05)
        #     for g in guesses:
        #         alt, az = self.moon_altaz(mjd=mjd+g)
        #         print(g, alt)
        #     raise Exception()

        if rise:
            next_brightness = np.float64(self.moon_illumination(mjd=event))
        else:
            next_brightness = np.float64(0)

        return np.float64(event), next_brightness

    def deltaV_sky_pos(self, mjd, targ_ra, targ_dec):
        """create inputs to KS91 from convenient params, return deltaV array

        Parameters:
        ----------
        mjd: np.float64
            decimal mjd to compute delta for

        targ_ra: np.float64 (or array)
            target ra

        targ_dec: np.float64 (or array)
            target dec

        Returns:
        -------
        deltaV : np.float64 (or array)
            change in V mag at targ location due to moon
        """

        moon_pos = self.moon_radec(mjd=mjd)
        lunar_phase = self.moon_illumination(mjd=mjd)

        alpha = 180.*np.arccos(2*lunar_phase-1)/3.1415

        moon_targ_dist = pyasl.getAngDist(moon_pos[0], moon_pos[1], targ_ra, targ_dec)

        malt, maz = self.radec2altaz(mjd=mjd, ra=moon_pos[0], dec=moon_pos[1])

        talt, taz = self.radec2altaz(mjd=mjd, ra=targ_ra, dec=targ_dec)

        zee = 90 - talt

        zee_m = 90 - malt

        # print(moon_pos)

        deltaV = KS91_deltaV(alpha, moon_targ_dist, zee, zee_m)

        # if len(moon_targ_dist) > 1:
        #     for m, d in zip(moon_targ_dist, deltaV):
        #         print(mjd, f"{float(lunar_phase):.2f} {float(alpha):3.1f} {float(m):.1f} {d}")
        #     for tt, z in zip(talt, zee):
        #         print(f"{float(tt):.1f} *{float(z):.1f}* {float(malt):.1f} *{float(zee_m):.1f}*")
        # else:
        #     print(mjd, f"{float(lunar_phase):.2f} {float(alpha):3.1f} {float(moon_targ_dist):.1f} {deltaV}")

        return deltaV

    def next_change(self, mjd):
        """How long before sky brightness changes drastically?
           Keep from observing a field too long and stealing dark/bright time

        Parameters:
        ----------
        mjd: np.float64
            decimal mjd to compute next for

        Returns:
        -------
        change : np.float64 (or array)
            mjd of next change
        """

        e_twi = self.evening_twilight(int(np.round(mjd)), twilight=self.dark_twilight)

        if mjd < e_twi:
            brightness = self.skybrightness(mjd + 1 / 24)
            if brightness < 0.35:
                return e_twi, np.float64(brightness)

        change, next_brightness = self.moon_rise_set(mjd)

        twi = self.morning_twilight(int(np.round(mjd)), twilight=self.dark_twilight)

        if twi < change:
            change = np.float64(twi)
            next_brightness = np.float64(1)

        if twi < mjd:
            twi = self.morning_twilight(int(np.round(mjd)), twilight=self.bright_twilight)
            change = np.float64(twi)
            next_brightness = np.float64(1)

        return change, next_brightness


class Master(Observer):
    """Master class to interpret master schedule as an observer. Inherits from Observer.

    Parameters:
    ----------
    schedule : str
        schedule file to use; default $ROBOSCHEDULER_DIR/data/master_schedule.par

    Attributes:
    ----------
    start : np.int32
        MJD (days) of first night of survey
    end : np.int32
        MJD (days) of last night of survey
    mjds : ndarray of np.int32
        MJDs (days) when survey is potentially active
    events : ndarray of numpy.str_
        names of events of note
    event_dates : ndarray of numpy.str_
        list of dates in ISO format for events of note
    event_times : ndarray of numpy.str_
        list of times in ISO format for events of note
    event_mjd : ndarray of numpy.float64
        MJDs (days) of events of note
    """
    def __init__(self, schedule='normal', observatory='apo',
                 observatoryfile=None):
        """Create Master object for schedule"""
        super().__init__(observatory=observatory,
                         observatoryfile=observatoryfile)
        masterfile = 'master_schedule_{o}_{s}.par'.format(o=observatory,
                                                          s=schedule)
        prod_dir = os.path.abspath(__file__).split("/scheduler.py")[0]
        schedulefile = os.path.join(prod_dir,
                                    'etc', masterfile)
        self._schedulefile = schedulefile
        self.schedule = yanny.yanny(self._schedulefile)
        self._validate()
        self.event_dates = np.array([date.decode() for date
                                     in self.schedule['SCHEDULE']['date']])
        self.event_times = np.array([time.decode() for time
                                     in self.schedule['SCHEDULE']['time']])
        self.event_mjds = self._dateandtime2mjd()
        self.events = np.array([event.decode() for event
                                in self.schedule['SCHEDULE']['event']])
        self.start = self._start()
        self.end = self._end()
        self.mjds = self._mjds()
        self.dark_twilight = np.float32(self.schedule['dark_twilight'])
        self.bright_twilight = np.float32(self.schedule['bright_twilight'])
        return

    def _dateandtime2mjd(self):
        return(dateandtime2mjd(date=self.event_dates,
                               time=self.event_times,
                               to_tai=self.schedule['to_tai']))

    def _validate(self):
        # should make sure:
        #  one start (first event)
        #  one end (last event)
        #  start MJD is a daytime time
        #  START_SURVEY is "on"
        #  END_SURVEY is "off"
        return

    def on(self, mjd=None):
        if(mjd < self.event_mjds[0]):
            return('off', self.event_mjds[0])
        if(mjd >= self.event_mjds[-1]):
            return('off', mjd + 1.)
        # Assumes there is only one
        indx = np.where((mjd >= self.event_mjds[0:-1]) &
                        (mjd < self.event_mjds[1:]))[0][0]
        return(self.schedule[self.events[indx]],
               self.event_mjds[indx + 1])

    def end_mjd(self):
        """Return end MJD

        Returns:

        end_mjd : np.float64
            MJD of last event (end of survey)
        """
        return(self.event_mjds[-1])

    def _start(self):
        # Assumes there is only one
        indx = np.where(self.events == 'START_SURVEY')[0][0]
        # Assumes START_SURVEY turns state on
        return(np.int32(np.floor(self.event_mjds[indx])))

    def _end(self):
        # Assumes there is only one
        indx = np.where(self.events == 'END_SURVEY')[0][0]
        # Assumes END_SURVEY turns state off
        return(np.int32(np.ceil(self.event_mjds[indx])))

    def _mjds(self):
        nmjd = self.end - self.start + 1
        mjds = self.start + np.arange(nmjd, dtype=np.int32)
        keep = np.zeros(nmjd, dtype=np.int32)
        for indx in np.arange(len(self.events) - 1):
            this_event = self.events[indx]
            if(self.schedule[this_event] == 'on'):
                keep_start = np.int32(np.floor(self.event_mjds[indx]))
                keep_end = np.int32(np.ceil(self.event_mjds[indx + 1]))
                ikeep = np.where((mjds >= keep_start) &
                                 (mjds <= keep_end))[0]
                keep[ikeep] = 1
        ikeep = np.where(keep)[0]
        return(mjds[ikeep])


class Scheduler(Master):
    """Scheduler class to keep track of fields and schedule the optimal field at
    a given time. Inherits from Master.

    Parameters:
    ----------

    airmass_limit : float, np.float32
        airmass limit for observations

    observatory : str
        Name of observatory to use (must be in observatory file)
        (default 'apo')

    observatoryfile : str
        Name of Yanny-format observatory file to read
        (default $ROBOSCHEDULER_DIR/data/observatories.par)

    Attributes:
    ----------

    airmass_limit : float, np.float32
        airmass limit for observations

    fields : Fields object
        object for fields

    observations : Observations object
        object accessing list of observations

    Comments:
    --------

    Scheduling proceeds conceptually as follows
         - fields are limited to set that are conceivably observable
         - A strategy to optimize completion

    """
    def __init__(self, airmass_limit=2.,
                 schedule='normal', observatory='apo', observatoryfile=None,
                 exp_time=None):
        """Return Scheduler object
        """
        super().__init__(schedule=schedule, observatory=observatory,
                         observatoryfile=observatoryfile)
        self.airmass_limit = airmass_limit
        if exp_time is None:
            self.exp_time = 18 / 60 / 24
        else:
            self.exp_time = exp_time
        return

    def initdb(self, designbase='plan-0', fromFits=True):
        """Initialize Scheduler fields and observation lists.
        Required before fields can be scheduled.

        Parameters:
        ----------

        designbase : str
            the name of the robostrategy version to load.

        fromFits : boolean
            For simulations we want to use a fits file, which should
            exist in "$OBSERVING_PLAN_DIR" and be named using 'designbase'
        """
        self.cadencelist = roboscheduler.cadence.CadenceList()
        self.fields = roboscheduler.fields.Fields(plan=designbase,
                                                  observatory=self.observatory)
        if fromFits:
            filebase = os.path.join(os.getenv('OBSERVING_PLAN_DIR'),
                                    designbase)
            # base = os.getenv('OBSERVING_PLAN_DIR')
            cadence_file = filebase + "/" + "rsCadences" + "-" + designbase + "-"\
                           + self.observatory + ".fits"
            fields_file = filebase + "/" + "rsAllocation" + "-" + designbase + "-"\
                           + self.observatory + ".fits"

            self.cadencelist.fromfits(filename=cadence_file)
            self.fields.fromfits(filename=fields_file)
        else:
            self.cadencelist.fromdb()
            self.fields.fromdb()

        self.observations = roboscheduler.observations.Observations(observatory=self.observatory)
        return

    def observable(self, mjd=None,  maxExp=None, check_skybrightness=True,
                   check_cadence=True, ignore=[]):
        """Return array of fields observable

        Parameters:
        ----------

        mjd : np.float64
            current MJD

        maxExp : integer
            the maximum number of exposures to allow, useful towards

        check_skybrightness : boolean
            passed to cadence check, rarely used

        check_cadence : boolean
            if all else fails, just see if a field is up, don't worry
            about cadence, rarely used

        ignore : list
            a list of fields to mark unobservable. Mostly used while
            planning a night to avoid rescheduling the same field that
            hasn't been marked done in the database.

        """

        (alt, az) = self.radec2altaz(mjd=mjd, ra=self.fields.racen,
                                     dec=self.fields.deccen)
        skybrightness = self.skybrightness(mjd)
        # valid cadence checks against "none" cadence issue
        observable = (alt > 0.) & self.fields.validCadence & self.fields.notDone
        nexp = np.ones(len(observable), dtype=int)
        delta_priority = np.zeros(len(observable), dtype=np.float64)

        deltav = self.deltaV_sky_pos(mjd, self.fields.racen, self.fields.deccen)
        airmass = self.alt2airmass(alt)
        moon_dist = self.moon_dist(mjd=mjd, ra=self.fields.racen,
                                   dec=self.fields.deccen)

        # whereRM = np.where(["bhm_rm" in c for c in self.fields.cadence])[0]

        next_change, next_brightness = self.next_change(mjd)

        nexp_change = int((next_change - mjd) / self.exp_time)

        if next_brightness <= 0.35:
            if skybrightness <= 0.35:
                # it's still dark
                nexp_change = 1e3
        else:
            # next is bright
            if skybrightness > 0.35:
                # it's bright now, don't change
                nexp_change = 1e3
            elif nexp_change == 0:
                # don't waste time, make it bright
                nexp_change = 1
                skybrightness = next_brightness

        # print(f"{float(mjd):.3f} {float(next_change):.3f} {float(next_brightness):.2f} {nexp_change}", maxExp)

        # indxs = np.where(self.fields.nextmjd > mjd)[0]
        # observable[indxs] = False
        indxs = np.where(observable)[0]
        for indx in indxs:
            # if(observable[indx]):
            if int(self.fields.field_id[indx]) in ignore:
                observable[indx] = False
                continue
            elif self.fields.flag[indx] == -1:
                observable[indx] = False
                continue
            cadence = self.cadencelist.cadences[self.fields.cadence[indx]]

            mjd_past = self.fields.hist[self.fields.field_id[indx]]
            # epoch_idx is the *index* of the *next* epoch
            # for 0 indexed arrays, this equivalent to
            # "how many epochs have I done previously"
            tol = cadence.max_length[int(self.fields.epoch_idx[indx])]
            epoch_idx, mjd_prev = epochs_completed(mjd_past, tolerance=tol)
            # how many exp/"designs" since start of last epoch?
            exp_epoch = np.sum(np.greater(mjd_past, mjd_prev))
            if exp_epoch < cadence.nexp[epoch_idx] and exp_epoch != 0:
                nexp[indx] = cadence.nexp[epoch_idx] - exp_epoch
                epoch_idx -= 1
                partial_epoch = True
            else:
                nexp[indx] = cadence.nexp[epoch_idx]
                partial_epoch = False
            if epoch_idx >= cadence.nepochs:
                print("DONE ", epoch_idx, cadence.nepochs, int(self.fields.field_id[indx]))
                observable[indx] = False
                continue
            if nexp[indx] > nexp_change and self.fields.flag[indx] != 1:
                # change between bright/dark, field doesn't fit
                observable[indx] = False
                continue
            observable[indx], delta_priority[indx] =\
                cadence.evaluate_next(epoch_idx=epoch_idx,
                                      partial_epoch=partial_epoch,
                                      mjd_past=mjd_prev,
                                      mjd_next=mjd,
                                      skybrightness_next=skybrightness,
                                      moon_dist=moon_dist[indx],
                                      deltaV=deltav[indx],
                                      airmass=airmass[indx])
            # if nexp[indx] > 6 and observable[indx]:
            #     delta_priority[indx] += 1e6
            delta_priority[indx] += nexp[indx] * 10
            if nexp[indx] > maxExp:
                observable[indx] = False
            if self.fields.flag[indx] == 1:
                # flagged as top priority
                # if we're in this loop, it's above the horizon
                # moon may not be ok, caveat emptor
                # so override cadence eligibility and bump priority
                observable[indx] = True
                delta_priority[indx] += 1e6
        # print(f"{float(mjd):.3f} {float(next_change):.3f} {float(next_brightness):.2f} {nexp_change}", maxExp)
        iobservable = np.where(observable)[0]

        return iobservable, nexp[iobservable], delta_priority[iobservable]

    def prioritize(self, mjd=None, iobservable=None, nexp=None,
                   delta_priority=None):
        """Prioritize fields according to the Robostrategy LST plan
        and accounting for priority adjustments from the cadence check (i.e.
        whether a field is inside an incomplete epoch or else how long
        before it fails cadence requirements.)

        Parameters:
        ----------

        mjd : float
            mjd to check, for find current lst
        iobservable : list of integer
            index into Scheduler.fields of observable fields
        nexp : ndarray  of np.int32
            array of number of exp for each field
        delta_priority : ndarray of np.float64
            change to the default priority, from cadence check.
            This accounts for time left in cadence or whether
            the cadence is partially complete.

        Returns:
        -------

        priority : ndarray of np.float64
            calculated priority for each field
        """

        priority = np.ones(len(iobservable))*100
        # priority = self.fields.basePriority[fieldid]
        priority += delta_priority

        priority += np.power(2, nexp) * 5  # 1280 for RM

        lst = self.lst(mjd)

        lstHrs = lst/15

        # lstDiffs = lstDiff(self.fields.lstPlan[fieldid], np.ones(len(fieldid))*lstHrs)

        lstDiffs = self.fields.lstWeight(lstHrs, iobservable)

        assert len(lstDiffs) == len(iobservable), "lst weight going poorly"
        # assert 0 not in delta_remaining, "some delta remaining not set properly!"

        ha = self.ralst2ha(ra=self.fields.racen[iobservable], lst=lst)
        dec = self.fields.deccen[iobservable]

        # gaussian weight, mean already 0, use 1 hr  std
        priority += 50 * np.exp(-(lstDiffs)**2 / (2 * 0.5**2))
        # gaussian weight, mean already 0, use 1 hr = 15 deg std
        priority += 50 * np.exp(-(ha)**2 / (2 * 15**2))
        # gaussian weight, mean = obs lat, use 20 deg std
        priority -= 20 * np.exp(-(dec - self.latitude)**2 / (2 * 20**2))

        return priority

    def pick(self, priority=None, fieldid=None, nexp=None):
        """Choose the next field

        Parameters:
        ----------

        priority : ndarray of np.float64
            calculated priority for each field
        fieldid : ndarray of np.int32
            field ids corresponding to prioritiea
        nexp : ndarray  of np.int32
            array of number of exp for each field

        Returns:
        -------

        pick_fieldid : np.int32
            the field id corresponding to the highest priority
        pick_exp : integer
            the number of exposures needed in the chosen field

        """
        assert len(priority) == len(fieldid) and len(priority) == len(nexp), \
            "inputs must be same size!"
        ipick = np.argmax(priority)
        pick_fieldid = fieldid[ipick]
        pick_exp = nexp[ipick]

        # if self.check:
        #     self.chosen[-1] = pick_fieldid
        #     self.chosen_pri[-1] = float(priority[ipick])
        #     print("CHOICE", pick_fieldid, np.max(priority), ipick, pick_exp, "\n")

        return(pick_fieldid, pick_exp)

    def designsNext(self, fieldid):
        """Figure out next designs on a field, i.e. which exposure are we on?

        Parameters:
        ----------

        fieldid : ndarray of np.int32
            field ids corresponding to prioritiea

        Returns:
        -------

        designs : list of integer
            a list of exposure numbers up next for "fieldid", corresponding 
            to some designs.

        """
        fieldidx = self.fields.getidx(fieldid)
        mjd_past = self.fields.hist[fieldid]
        cadence = self.cadencelist.cadences[self.fields.cadence[fieldidx]]
        epoch_idx, mjd_prev = epochs_completed(mjd_past, tolerance=240)
        nexp_next = cadence.nexp[epoch_idx]
        # how many completed designs since start of epoch
        # this assumes designs are marked incomplete when epoch_max_length
        # is hit
        exp_in_epoch = np.sum(np.greater(mjd_past, mjd_prev))

        exp_to_do = nexp_next - exp_in_epoch

        designs = list(len(mjd_past) + np.arange(exp_to_do))

        return designs

    def nextfield(self, mjd=None, maxExp=None, returnAll=False, live=False,
                  ignore=[]):
        """Picks the next field to observe

        Parameters:
        ----------

        mjd : np.float64
            Current MJD (days)
        maxExp : int
            maximum number of full exposures before next event
        returnAll : boolean
            Return all the fields? For choosing backups
        live : boolean
            Are we live in Kronos? otherwise sim behavior slightly different
        ignore : list
            Fields to ignore temporarily, i.e. because they are already
            scheduled but have no exposures in opsdb

        Returns:
        --------

        fieldid : np.int32, integer
            ID of field to observe
        designs : list of integer
            a list of exposure numbers up next for "fieldid", corresponding
            to some designs.

            OR if returnAll

        sorted_fields : list of integer
            a list of fields sorted by priority
        sorted_exp : list of integer
            a list of nexp sorted by priority
        """

        iobservable, nexp, delta_priority = self.observable(mjd=mjd, maxExp=maxExp,
                                                            ignore=ignore)
        if(len(iobservable) == 0) and live:
            # in a sim we'll take the dead time, should write something to track this better
            # print("Nothing observable")
            iobservable, nexp, delta_priority = self.observable(mjd=mjd, maxExp=maxExp,
                                                                ignore=ignore)
        if len(iobservable) == 0:
            # print("!! nothing to observe; {} exp left in the night".format(maxExp))
            if returnAll:
                return None, -1
            return None, -1

        priority = self.prioritize(iobservable=iobservable, mjd=mjd,
                                   nexp=nexp, delta_priority=delta_priority)

        if returnAll:
            # inverse priority, highest first
            sorted_priority = np.argsort(priority)[::-1]
            # for i in sorted_priority:
            #     field_idx = iobservable[i]
            #     print(i, priority[i], nexp[i],
            #           self.fields.field_id[field_idx],
            #           self.fields.cadence[field_idx])
            sorted_idx = [iobservable[i] for i in sorted_priority]
            sorted_fields = [self.fields.field_id[i] for i in sorted_idx]
            sorted_exp = [nexp[i] for i in sorted_priority]

            if not live:
                # for observesim for now I'm sorry!
                return np.array(sorted_idx), np.array(sorted_exp)

            return sorted_fields, sorted_exp

        # considered = False
        # print(observable_fieldid)
        # print(priority, self.fields.cadence[observable_fieldid])
        # for p, c, i in zip(priority, np.array(self.fields.cadence)[observable_fieldid], observable_fieldid):
        #     if "RM" in c.upper():
        #         print(c, i, p, np.max(priority))
        #         considered = True

        observable_fieldid = self.fields.field_id[iobservable]

        fieldid, next_exp = self.pick(priority=priority,
                                      fieldid=observable_fieldid,
                                      nexp=nexp)

        if not live:
            # just a sim return number of designs
            return(fieldid, next_exp)

        # its live, return list of exp indices
        designs = self.designsNext(fieldid)

        return fieldid, designs

    def update(self, fieldid=None, result=None, finish=False):
        """Update Scheduler.observations with result of observations, used for sims

        Parameters:
        -----------

        fieldid : np.int32, int
            ID of field

        result : ndarray
            One element, contains 'mjd', 'duration', 'sn2'

        """

        fieldidx = int(np.where(self.fields.field_id == fieldid)[0])

        racen = self.fields.racen[fieldidx]
        deccen = self.fields.deccen[fieldidx]
        cadence = self.fields.cadence[fieldidx]
        nfilled = self.fields.nfilled[fieldidx]
        nexp_cumul = len(self.fields.observations[fieldidx]) + 1

        (alt, az) = self.radec2altaz(mjd=result['mjd'],
                                     ra=racen,
                                     dec=deccen)
        airmass = self.alt2airmass(alt)
        skybrightness = self.skybrightness(result['mjd'])
        lst = self.lst(result['mjd'])
        iobs = self.observations.add(fieldid=fieldid,
                                     mjd=result['mjd'],
                                     duration=result['duration'],
                                     apgSN2=result['apgSN2'],
                                     rSN2=result['rSN2'],
                                     bSN2=result['bSN2'],
                                     skybrightness=skybrightness,
                                     airmass=airmass,
                                     lst=lst,
                                     racen=racen,
                                     deccen=deccen,
                                     cadence=cadence,
                                     nfilled=nfilled,
                                     nexp_cumul=nexp_cumul)

        # iobservations = self.fields.observations[fieldidx]
        # mjd_past = self.observations.mjd[iobservations]
        # epoch_idx is the *index* of the *next* epoch
        # for 0 indexed arrays, this equivalent to
        # "how many epochs have I done previously"
        # epoch_idx, mjd_prev = epochs_completed(mjd_past, tolerance=45)
        if finish:
            # if fieldidx == 4906:
            #     print(self.fields.hist[fieldid])
            #     print("UPDATE", float(result['mjd']))
            self.fields.completeDesign(fieldidx, float(result['mjd']), lst, iobs)
        # self.fields.add_observations(result['mjd'], fieldidx, iobs, lst,
        #                              epoch_idx)
        return


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
