# Some variables that will soon be moved to a separate file
TARGET_ELEVATION_MIN = 20 # this elevation is the physical minimum, below this the ADC does not work
TARGET_ELEVATION_HIGH_MIN = 45 # this elevation is the preferred one for stars that will be high in the sky
TARGET_ELEVATION_MAX = 85
TARGET_EXPOSURE_TIME_MAX =  2* 60 * 60 # 2 hour
TARGET_MOON_DIST_MIN = 15
TARGET_MOON_DIST_MAX = 25

# Maximum single exposure time in seconds
MAX_EXPTIME = 1200.
MIN_EXPTIME = 300.
MIN_TOTOBS = 300.
READOUT = 40.

MAX_I2 = 60000
MIN_I2 = 500

MAX_EXPMETER = 2e9 # above this saturates the CCD

MAX_NSHOTS = 10 # because people are stupid

# A few constants to make accessing the star table more readable
DS_RA     = 0
DS_DEC    = 1
DS_PMRA   = 2
DS_PMDEC  = 3
DS_VMAG   = 4
DS_EXPT   = 5
DS_COUNTS = 6
DS_APFPRI = 7
DS_CAD    = 8
DS_NSHOTS = 9
DS_LAST   = 10
DS_BV     = 11
DS_ERR    = 12
DS_UTH    = 13
DS_UTM    = 14
DS_DUR    = 15
DS_MIN    = 16
DS_MAX    = 17
DS_NOB    = 18
DS_TOT    = 19
DS_I2CNTS = 20
DS_SHEETN = 21

SLOWDOWN_MIN = 0.6
SLOWDOWN_THRESH = 2.5
SLOWDOWN_MAX = 5.0
SLOWDOWN_VMAG_LIM = 9.0 
SEEING_THRESH = 30.0

PRI_DELTA = 20

DS_FT_NAMES = 0
DS_FT_FRAC = 1
DS_FT_TOT = 2
DS_FT_CUR = 3

EXP_LIM = 3e9
