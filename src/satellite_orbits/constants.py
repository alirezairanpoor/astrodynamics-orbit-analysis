# Earth constants used throughout the project.
# Units: km, s

MU_EARTH = 398600.4405  # km^3 / s^2
R_EARTH = 6378.137      # km
J2_EARTH = 1.08263e-3
# WGS84 reference ellipsoid
FLATTENING_EARTH = 1 / 298.257223563

# Speed of light in vacuum.
# Exact by the SI definition of the metre.
C_LIGHT = 299792.458  # km / s

# Earth rotation rate with respect to the inertial frame.
# This is the time derivative of Greenwich Mean Sidereal Time,
# not 2*pi/86400, because the sidereal day is shorter than the
# solar day by about 236 s.
OMEGA_EARTH = 7.2921158553e-5  # rad / s

# Time conversions.
SECONDS_PER_DAY = 86400.0  # s
# GPS carrier frequencies.
GPS_L1_FREQUENCY = 1575.42e6  # Hz
GPS_L2_FREQUENCY = 1227.60e6  # Hz
