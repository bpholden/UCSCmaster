from __future__ import print_function

from datetime import datetime, timedelta
import os
import sys
import time
import pytz

import numpy as np
import astropy
import astropy.coordinates
import astropy.units
import astroplan

import SchedulerConsts


def makeAPFObs():
    # Generate a astropy.coordinate observer for the APF

    apf_lat = astropy.coordinates.Latitude((37,20,33.1),unit=astropy.units.deg)
    apf_long = astropy.coordinates.Longitude((-121,38,17.7),wrap_angle=astropy.units.deg*180,unit=astropy.units.deg)
    apf_height = 1274 * astropy.units.meter
    apf_loc = astropy.coordinates.EarthLocation.from_geodetic(apf_long,apf_lat,apf_height)

    apf_obs = astroplan.Observer(name='APF Telescope',
               location=apf_loc,
               pressure=870 * astropy.units.hPa,
               relative_humidity=0.3,
               temperature=20 * astropy.units.deg_C,
               timezone=pytz.timezone('US/Pacific'),
               description="APF Telescope on Mount Hamilton, California")
    
    return apf_obs



def visible(observer, cdate, stars, obs_lens, pref_min_el=SchedulerConsts.TARGET_ELEVATION_HIGH_MIN, min_el=SchedulerConsts.TARGET_ELEVATION_MIN,
                   max_el=SchedulerConsts.TARGET_ELEVATION_MAX):
    """ Args:
            stars: A list of pyephem bodies to evaluate visibility of
            observer: A pyephem observer to use a the visibility reference
            obs_len: A list of observation lengths ( Seconds ). This is the time frame for which visibility is checked

            Optional:
            pref_min_el: Preferred minimum body elevation to be visible ( degrees )            
            min_el: The minimum body elevation to be visible ( degrees ) - only use this if star never goes above preferred limit
            max_el: The maximum body elevation to be visible ( degrees )
        Returns:
            Boolean list representing if body[i] is visible

        Notes: Uses the observer's current date and location
    """

    ret = []
    fin_elevations = []
    start_elevations = []

    
    # Now loop over each body to check visibility
    for star, obslen in zip(stars, obs_lens):

        constraints = [astroplan.AltitudeConstraint(min_el*astropy.units.deg, max_el*astropy.units.deg)]        

        if obs_len> 0:
            findate = cdate + timedelta(seconds=obs_len)
        else:
            findate = cdate + timedelta(seconds=1)
            
        time_range = astropy.time.Time([cdate,findate])

        altaz = apf_obs.altaz(cdate,target=star)
        cur_el = altaz.alt.value
        start_elevations.append(cur_el)

        fin_altaz = apf_obs.altaz(findate,target=star)
        fin_el = fin_altaz.alt.value
        fin_elevations.append(fin_el)
        
        rv = astroplan.is_always_observable(constraints, apf_obs, star, time_range=time_range)
        if len(rv) == 1:
            ret.append(rv[0])
        else:
            ret.append(False)
        

    return ret, np.array(start_elevations), np.array(fin_elevations)



def visibleSE(observer, cdate, stars, obs_len, pref_min_el=SchedulerConsts.TARGET_ELEVATION_HIGH_MIN, min_el=SchedulerConsts.TARGET_ELEVATION_MIN,
                   max_el=SchedulerConsts.TARGET_ELEVATION_MAX,shiftwest=False):
    """ Args:
            stars: A list of pyephem bodies to evaluate visibility of
            observer: A pyephem observer to use a the visibility reference
            obs_len: A list of observation lengths ( Seconds ). This is the time frame for which visibility is checked
            pref_min_el: Preferred minimum body elevation to be visible ( degrees )            
            min_el: The minimum body elevation to be visible ( degrees ) - only use this if star never goes above preferred limit
            max_el: The maximum body elevation to be visible ( degrees )
        Returns:
            Boolean list representing if body[i] is visible

        Notes: Uses the observer's current date and location
    """
    # Store the previous observer horizon and date since we change these
    prev_horizon = observer.horizon
    cdate = observer.date
    ret = []
    fin_elevations = []
    start_elevations = []
    scaled_elevations = []

    sun_alt_az = apf_obs.sun_altaz(cdate)

    sun_el = sun_alt_az.alt.value
    sun_az = sun_alt_az.az.value
    
    bottom_angle = SchedulerConsts.SUNEL_STARTLIM-15 # typically -24 degrees
    
    if sun_el > (bottom_angle) and sun_az > 180 and shiftwest:
        offset = 3*(sun_el - bottom_angle) # note, this is positive
        preferred_angle = (90 - offset)
    else:
        offset = 0.0
        preferred_angle = 90 

    
    # Now loop over each body to check visibility
    for star, obs_len in zip(stars, obs_lens):

        # Is the target visible now?
        constraints = [astroplan.AltitudeConstraint(min_el*astropy.units.deg, max_el*astropy.units.deg)]        

        if obs_len> 0:
            findate = cdate + timedelta(seconds=obs_len)
        else:
            findate = cdate + timedelta(seconds=1)
            
        time_range = astropy.time.Time([cdate,findate])

        altaz = apf_obs.altaz(cdate,target=star)
        cur_el = altaz.alt.value
        start_elevations.append(cur_el)

        fin_altaz = apf_obs.altaz(findate,target=star)
        fin_el = fin_altaz.alt.value
        fin_elevations.append(fin_el)
        
        rv = astroplan.is_always_observable(constraints, apf_obs, star, time_range=time_range)
        if len(rv) == 1:
            ret.append(rv[0])
        else:
            ret.append(False)

        diff = (star.dec.value - apf_obs.location.lat.value)
        transit_alt = 90.0 - diff
        se = 90.0 - (transit_alt - cur_el) 
        if offset > 0:
            if cur_az < 180:
                se -= offset
            else:
                se = 90 - np.abs(preferred_angle - se)
        scaled_elevations.append(se)
        
    return ret, np.array(start_elevations), np.array(fin_elevations), np.array(scaled_elevations)



if __name__ == '__main__':
    # Generate a pyephem observer for the APF
    apf_obs = makeAPFObs()
    cur_date = datetime.utcfromtimestamp(int(time.time()))

    star = astropy.coordinates.SkyCoord("1h44m4.083s","-15d56m14.93s")
    ret, se, fe = visible(apf_obs,cur_date,[star],[0.])
    print (ret, se, fe)
    ret, se, fe = visible(apf_obs,cur_date,[star],[400.])
    print (ret, se, fe)

    
    ret, se, fe, sce = visibleSE(apf_obs,cur_date,[star],[400.])
    print (ret, se, fe, sce)
