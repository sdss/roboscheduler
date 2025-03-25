#! /usr/bin/env python

from collections import defaultdict

import numpy as np
import scipy.optimize as optimize
from astropy.time import Time

from roboscheduler import scheduler


def _bright_dark_function(mjd=None, sched=None, switch=0.35):
        return sched.skybrightness(mjd) - switch


def summmerOrWinter_north(startTime):
    """Check whether we're between equinoxes, "winter"
       Expects datetime object to check
    """
    winter = startTime.month < 3
    if not winter and startTime.month == 3:
        winter = startTime.day <= 20
    fall = startTime.month >= 10
    if not fall and startTime.month == 9:
        fall = startTime.day >= 22
    return winter or fall


def nightSchedule(sched, night_start, night_end):
    fudge = 25 / 60 / 24
    bright_start = bool(sched.skybrightness(night_start + fudge) >= 0.35)
    bright_end = bool(sched.skybrightness(night_end - fudge) >= 0.35)
    dark_start = bool(sched.skybrightness(night_start + fudge) < 0.35)
    dark_end = bool(sched.skybrightness(night_end - fudge) < 0.35)

    # mjd_sched = dict()
    Bright_Start = 0
    Bright_End = 0
    Dark_Start = 0
    Dark_End = 0

    if bright_start and bright_end:
        Bright_Start = night_start
        Bright_End = night_end
    elif dark_start and dark_end:
        Dark_Start = night_start
        Dark_End = night_end
    elif dark_start and bright_end:
        split = optimize.brenth(_bright_dark_function,
                                night_start + fudge, night_end - fudge,
                                args=(sched, 0.35))
        Bright_Start = split
        Bright_End = night_end
        Dark_Start = night_start
        Dark_End = split
    elif bright_start and dark_end:
        split = optimize.brenth(_bright_dark_function,
                                night_start + fudge, night_end - fudge,
                                args=(sched, 0.35))
        Bright_Start = night_start
        Bright_End = split
        Dark_Start = split
        Dark_End = night_end
    
    dark_time = float(Dark_End - Dark_Start)
    bright_time = float(Bright_End - Bright_Start)
    return dark_time, bright_time


def mjd_dict():
    # we're abusing the ddict default_factory
    return {"bright": 0, "dark": 0, "twilight": 0}


def computeSched(loc):
    if loc == "lco":
        start = 60000
        end = start + 5.75*365.25
    else:
        start = 59640
        end = start + 5.0*365.25
    mjds = np.arange(start, end, 1)

    sched = scheduler.Scheduler(observatory=loc, schedule="v6")

    time_avail = defaultdict(mjd_dict)

    for m in mjds:
        onoff, nextchange_on = sched.on(mjd=m)

        if onoff != "on":
            print(m)
            continue

        mjd_evening_twilight = sched.evening_twilight(m, twilight=-15)
        mjd_morning_twilight = sched.morning_twilight(m, twilight=-15)
        
        atime = Time(m, format="mjd")
        atime.format = "datetime"

        winter = summmerOrWinter_north(atime.datetime)
        if loc == "lco":
            winter = not winter

        if winter:
            mjd_evening_twilight_bright = sched.evening_twilight(m, twilight=-12)
            mjd_morning_twilight_bright = sched.morning_twilight(m, twilight=-12)
        else:
            mjd_evening_twilight_bright = sched.evening_twilight(m, twilight=-8)
            mjd_morning_twilight_bright = sched.morning_twilight(m, twilight=-8)

        xtra_evening = mjd_evening_twilight - mjd_evening_twilight_bright
        xtra_morning = mjd_morning_twilight_bright - mjd_morning_twilight
        twilight = xtra_evening + xtra_morning
        
        dark_time, bright_time = nightSchedule(sched,
                                               mjd_evening_twilight,
                                               mjd_morning_twilight)
        
        time_avail[m]["bright"] = bright_time * 24
        time_avail[m]["dark"] = dark_time * 24
        time_avail[m]["twilight"] = twilight * 24

    with open(f"time_avail_{loc}.csv", "w") as of:
        cum_bright = 0
        cum_dark = 0
        cum_twilight = 0
        print("mjd, bright, cum_bright, dark, cum_dark, twilight, cum_twilight", file=of)
        for k,v in time_avail.items():
            cum_bright += v['bright']
            cum_dark += v['dark']
            cum_twilight += v['twilight']
            print(f"{k}, {v['bright']:.2f}, {cum_bright:.2f}, {v['dark']:.2f}, {cum_dark:.2f}, {v['twilight']:.2f}, {cum_twilight:.2f}", file=of)

if __name__ == "__main__":
    computeSched("lco")
    computeSched("apo")
