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


    ret = []
    
    obs_lens = np.asarray(obs_lens)
    
    constraints = [astroplan.AltitudeConstraint(min_el*astropy.units.deg, max_el*astropy.units.deg)]        

    sun_alt_az = observer.sun_altaz(cdate)

    sun_el = sun_alt_az.alt.value
    sun_az = sun_alt_az.az.value
    
    bottom_angle = SchedulerConsts.SUNEL_STARTLIM-15 # typically -24 degrees
    
    if sun_el > (bottom_angle) and sun_az > 180 and shiftwest:
        offset = 3*(sun_el - bottom_angle) # note, this is positive
        preferred_angle = (90 - offset)
    else:
        offset = 0.0
        preferred_angle = 90 

    obs_lens[obs_lens <=0] = 1.0

    start_dates = [ cdate for i in range(0,len(obs_lens))]
    fin_dates = start_dates + obs_lens / 86400.

    altaz = observer.altaz(start_dates,target=stars)
    start_elevations = altaz.alt.value
    cur_azs = altaz.az.value
    
    altaz = observer.altaz(fin_dates,target=stars)
    fin_elevations = altaz.alt.value

    diff = np.abs(stars.dec.value - observer.location.lat.value)
    transit_alts = 90.0 - diff
    scaled_elevations = 90.0 - (transit_alts - start_elevations)
    if offset > 0:
        scaled_elevations[cur_azs<180] -= offset
        scaled_elevations[cur_azs>180] = 90 - np.abs(preferred_angle - scaled_elevations[cur_azs>180])

    # Now loop over each body to check visibility
    for i in range(0,len(obs_lens)):

        # Is the target visible now?
        time_range = [start_dates[i],fin_dates[i]]
        rv = astroplan.is_always_observable(constraints, observer, stars, time_range=time_range,time_grid_resolution=10*astropy.units.second)
        
        if len(rv) > 0:
            ret.append(rv[i])
        else:
            ret.append(False)

        
    return np.asarray(ret), start_elevations, fin_elevations, scaled_elevations


if __name__ == '__main__':
    # Generate a pyephem observer for the APF

    apf_obs = makeAPFObs()
    cur_date = astropy.time.Time.now()
    
    star = astropy.coordinates.SkyCoord("1h44m4.083s","-15d56m14.93s")
    ret, se, fe, sce = visible(apf_obs,cur_date,star,np.asarray([0.]))
    print (ret, se, fe, sce)
    ret, se, fe, sce = visible(apf_obs,cur_date,star,np.asarray([400.]))
    print (ret, se, fe, sce)


