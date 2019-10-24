# Some variables that will soon be moved to a separate file
TARGET_ELEVATION_MIN = 20 # this elevation is the physical minimum, below this the ADC does not work
TARGET_ELEVATION_HIGH_MIN = 45 # this elevation is the preferred one for stars that will be high in the sky
TARGET_ELEVATION_MAX = 85
TARGET_EXPOSURE_TIME_MAX =  1* 60 * 60 # 1 hour
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

SLOWDOWN_MIN = 0.4
SLOWDOWN_MAX = 5.0

PRI_DELTA = 5

# A few constants to make accessing the star table more readable
TOO_RA     = 0
TOO_DEC    = 1
TOO_PMRA   = 2
TOO_PMDEC  = 3
TOO_VMAG   = 4
TOO_EXPT   = 5
TOO_COUNTS = 6
TOO_APFPRI = 7
TOO_CAD    = 8
TOO_NSHOTS = 9
TOO_LAST   = 10
TOO_UTH    = 11
TOO_UTM    = 12
TOO_DUR    = 13
TOO_NOB    = 14
TOO_TOT    = 15

IMMEDIATE = 3
NEXTEXP = 2
NEXTLINE = 1

PRI_IMMEDIATE = 300
PRI_NEXT_EXP = 200

