import numpy as np

def getRARad(hr, mn, sec):
    try:
        if hr < 0 or hr > 23:
            return None
        if mn < 0 or mn >= 59:
            return None
        if sec < 0 or sec >= 60:
            return None
        ra_hours = float(hr) + float(mn)/60. + float(sec)/3600.
        return ra_hours * 15 * np.pi/180.0
    except:
        return None

def getDECRad(deg, mn, sec, neg=False):
    try:
        deg = float(deg)
        mn = float(mn)
        sec = float(sec)
        if deg < -40 or deg > 90:
            return None
        if mn < 0 or mn >= 59:
            return None
        if sec < 0 or sec >= 60:
            return None
    except:
        return None
    if deg < 0:
        neg = True
        deg = abs(deg)       
    if  mn < 0:
        neg = True
        mn = abs(mn)
    if sec < 0:
        neg = True
        sec = abs(sec)
    x = deg + mn/60. + sec/3600.
    x = x * np.pi/180.
    if neg:
        return x*-1
    else:
        return x

def getCoordStr(floatval,isRA=False):

    neg = False
    nround = 2
    if isRA:
        floatval /= 15.
        nround = 3
    if floatval < 0:
        neg = True
    floatval = abs(floatval)
    deghrval = int(floatval)
    minval = (floatval % 1) * 60.0 
    secval = round( (minval % 1) *60.0, nround)

    if neg and deghrval != 0:
        ret = "-" + str(deghrval) + ' '
    else:
        ret = str(deghrval) + ' '
    if neg and deghrval == 0 and minval != 0:
        ret += "-" + str(int(minval)) + ' '
    else:
        ret += str(int(minval)) + ' '
    if neg and deghrval == 0 and minval == 0:
        ret += "-" + str(secval)
    else:
        ret += str(secval)
    return ret


def getLST(date, longitude):
    """Take a datetime and longitude and calculate the Local Sidereal Time."""
    # Assumes date is a datetime object, and that the longitude is formatted as in PyEphem 

    ll = [float(v) for v in longitude.split(':')]
    if ll[0] > 0:
        sign = 1
    else:
        sign = -1
    ut = date.hour + date.minute/60. + date.second/3600.
    lng = ll[0] + sign*ll[1]/60. + sign*ll[2]/3600.
    d  = ephem.julian_date() - 2451545.0
    lst = 100.46 + 0.985647 * d + lng + 15*ut
    return lst % 360.


        
def getElAz(ra, dec, lat, lng, time):
    """Given RA, DEC, Latitude, and a time, returns the corresponding elevation and azimuth angles
       Works with single values, or numpy arrays
       """
    lst = getLST(time, lng)
    ha = ((lst- np.degrees(ra)) % 360.) * np.pi/180.
    el = np.arcsin(np.sin(dec) * np.sin(lat) + \
                   np.cos(dec) * np.cos(lat) * np.cos(ha))
    az = np.arccos( (np.sin(dec) - np.sin(el)*np.sin(lat)) / \
                         (np.cos(el) * np.cos(lat)))
    return (np.degrees(el), np.degrees(az))

