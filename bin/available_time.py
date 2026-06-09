#! /usr/bin/env python

import argparse
from collections import defaultdict

import numpy as np
import scipy.optimize as optimize
from astropy.time import Time

from roboscheduler import scheduler

def sched_engineering(mjds, sched, apo=True):
    illum = sched.moon_illumination(mjds)

    full = np.where(illum > 0.98)[0]

    # days = list()
    nights = list()
    last_night = 0
    for i in full:
        m = int(mjds[i])
        if m - last_night < 5:
            continue
        aTime = Time(m-1, format="mjd").datetime
        if apo:
            # apo doesn't want weekends
            if aTime.isoweekday() == 6:
                nights.append(m-1)
                # days.append(m-2)
                last_night = m
            elif aTime.isoweekday() == 7:
                nights.append(m+1)
                # days.append(m+2)
                last_night = m
            elif aTime.isoweekday() == 1:
                nights.append(m)
                # days.append(m+1)
                last_night = m
            else:
                nights.append(m)
                # days.append(m-1)
                last_night = m
        else:
            # LCO doesn't want it on Mon->Tue
            if aTime.isoweekday() == 1:
                nights.append(m+1)
                last_night = m+1
            else:
                nights.append(m)
                last_night = m

    return nights


def plannedSkips(start, stop, loc="apo"):
    sched = scheduler.Observer(observatory=loc)

    apo = loc == "apo"

    mjds = np.arange(start, stop, 1)

    date = Time(mjds[0], format="mjd").datetime

    skipped_mjds = list()

    shutdown_duration = 7 * 6  # 6 weeks

    eng = sched_engineering(mjds, sched, apo)
    delay = 0
    for m in mjds[1:]:
        m = int(m)
        if m < delay:
            continue
        if m in eng:
            if apo:
                dur = 1
            else:
                dur = 2
            for i in range(dur):
                skipped_mjds.append(m)
            delay = m + dur
        elif apo:
            date = Time(m-1, format="mjd").datetime
            if date.month == 7:
                if date.isoweekday() == 1:
                    for i in range(shutdown_duration):
                        skipped_mjds.append(m+i)
                    delay = m + shutdown_duration
        else:
            # LCO
            date = Time(m-1, format="mjd").datetime
            if date.month == 12:
                if date.day == 24:
                    if (date.year - 2000) % 3 == 0:
                        # realuminization
                        dur = 10
                    else:
                        dur = 2
                    for i in range(dur):
                        skipped_mjds.append(m+i)
                    delay = m + dur

    return skipped_mjds


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
    split = None

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
    return dark_time, bright_time, split


def mjd_dict():
    # we're abusing the ddict default_factory
    return {"bright": 0, "dark": 0, "twilight": 0, 
            "lst_start": -1, "lst_end": -1, "lst_split": -1}


def computeSched(loc=None, start=None, end=None):
    mjds = np.arange(start, end, 1)

    sched = scheduler.Scheduler(observatory=loc, schedule="v6")

    time_avail = defaultdict(mjd_dict)

    skipped_mjds = plannedSkips(start, end, loc=loc)

    for m in mjds:
        if m in skipped_mjds:
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
        
        dark_time, bright_time, split = nightSchedule(sched,
                                                      mjd_evening_twilight,
                                                      mjd_morning_twilight)
        
        time_avail[m]["bright"] = bright_time * 24
        time_avail[m]["dark"] = dark_time * 24
        time_avail[m]["twilight"] = twilight * 24
        time_avail[m]["lst_start"] = float(sched.lst(mjd_evening_twilight))
        time_avail[m]["lst_end"] = float(sched.lst(mjd_morning_twilight))
        if split is not None:
            time_avail[m]["lst_split"] = float(sched.lst(split))

    with open(f"time_avail_{loc}.csv", "w") as of:
        cum_bright = 0
        cum_dark = 0
        cum_twilight = 0
        print("mjd, bright, slots_bright, dark, slots_dark, twilight, slots_twilight, lst_start, lst_end, lst_split", file=of)
        for k,v in time_avail.items():
            cum_bright += v['bright']
            cum_dark += v['dark']
            cum_twilight += v['twilight']
            print((
                  f"{k}, {v['bright']:.2f}, {int(v['bright'] * 3)}, "
                  f"{v['dark']:.2f}, {int(v['dark'] * 3)}, "
                  f"{v['twilight']:.2f}, {int(v['twilight'] * 3)}, "
                  f"{v['lst_start']:.1f}, {v['lst_end']:.1f}, {v['lst_split']:.1f}"
                  ) , file=of)

if __name__ == "__main__":
    usage = "make_figs"
    description = "Post processing for observesim, make figs and webpage"
    parser = argparse.ArgumentParser(description=description, usage=usage)
    parser.add_argument("-l", "--location", dest="location", type=str,
                        required=False, help="observatory location",
                        default="lco")
    parser.add_argument("-s", "--start", dest="start", type=str,
                        required=False, help="start MJD", default=None)
    args = parser.parse_args()
    location = args.location
    start = args.start

    if start is None:
        if location == "apo":
            start = 61406
        else:
            start = 61679
        end = start + 5*365.25

    computeSched(loc=location, start=start, end=end)
