# Exercise 2.6 — Initial Orbit Determination

## Objective

Determine the Keplerian orbit from two ground-station observations containing:

- Azimuth
- Elevation
- Range
- Observation time

The observations are first converted to inertial position vectors and then used in Gauss's two-position preliminary orbit determination method.

---

## Given Data

Ground-station ECEF coordinates:

X = 1344.143 km  
Y = 6068.601 km  
Z = 1429.311 km

### Observation 1

Date      : 1999-04-02  
Time      : 00:30:00  
Azimuth   : 132.67 deg  
Elevation : 32.44 deg  
Range     : 16945.450 km

### Observation 2

Date      : 1999-04-02  
Time      : 03:00:00  
Azimuth   : 123.08 deg  
Elevation : 50.06 deg  
Range     : 37350.340 km

The observation interval is:

delta_t = 2.5 h = 9000 s

---

## Solution Procedure

### 1. Ground-Station Coordinates

The given ECEF position of the ground station is converted to geodetic longitude and latitude.

r_station_ecef
        ↓
ecef_to_geodetic()
        ↓
longitude / latitude / height

The longitude and latitude define the local East-North-Zenith frame of the station.

---

### 2. Observation to Local ENZ Vector

Each azimuth, elevation, and range observation is converted into a local station-to-satellite vector:

Azimuth / Elevation / Range
            ↓
az_el_range_to_enz()
            ↓
          s_enz

with:

s_east   = range × cos(elevation) × sin(azimuth)

s_north  = range × cos(elevation) × cos(azimuth)

s_zenith = range × sin(elevation)

---

### 3. ENZ to ECEF

The local station-to-satellite vector is rotated into the Earth-fixed frame:

s_enz
  ↓
enz_to_ecef()
  ↓
s_ecef

The satellite position relative to the center of the Earth is then:

r_satellite_ecef
=
r_station_ecef + s_ecef

---

### 4. Observation Time and Earth Rotation

Each calendar observation time is converted to Modified Julian Date:

calendar date
    ↓
calendar_to_mjd()
    ↓
   MJD

The Modified Julian Date is then used to compute the Greenwich Hour Angle:

MJD
 ↓
mjd_to_greenwich_hour_angle()
 ↓
Theta

---

### 5. ECEF to ECI

The Earth-fixed satellite position is transformed into the inertial frame:

r_ecef
  ↓
ecef_to_eci()
  ↓
r_eci

The resulting inertial position vectors are approximately:

r_eci_1 =
[ 11959.978
 -16289.478
  -5963.827 ] km

r_eci_2 =
[ 39863.390
 -13730.547
  -4862.350 ] km

---

## Gauss Two-Position Orbit Determination

The two inertial position vectors and their time separation are passed to the general Gauss orbit determination routine:

r_eci_1
r_eci_2
delta_t
   ↓
gauss_two_position_elements()
   ↓
a / e / i / RAAN / AOP / M

The Gauss method used here is implemented independently of this exercise and can therefore be reused with other pairs of position vectors.

---

### 6. Orbital Plane

The two position vectors determine the orbital plane.

The first in-plane unit vector is:

e_a = r_a / |r_a|

A second perpendicular in-plane direction is constructed from the second position.

Their cross product gives the unit vector normal to the orbital plane:

W = e_a × e_0

From W, the following elements are obtained:

i

RAAN

The argument of latitude at the first observation is also determined from the geometry.

---

### 7. Sector-to-Triangle Ratio

The geometry alone does not determine the complete ellipse.

The observation time interval is therefore used through Gauss's sector-to-triangle ratio:

r_a
r_b
delta_t
  ↓
gauss_sector_triangle_ratio()
  ↓
eta

The ratio eta is solved iteratively using the secant method.

---

### 8. Semi-Latus Rectum

With the sector-to-triangle ratio, the semi-latus rectum is obtained from:

p = (2 × triangle_area × eta / tau)^2

where:

tau = sqrt(mu) × delta_t

---

### 9. Remaining Orbital Elements

Using p and the geometry of the two position vectors, the following quantities are determined:

e

TA

The argument of perigee follows from:

AOP = argument_of_latitude - TA

The semi-major axis is:

a = p / (1 - e^2)

The true anomaly is converted to eccentric anomaly:

TA → E

and finally the mean anomaly at the first observation is:

M = E - e × sin(E)

---

## Final Result

The resulting Keplerian elements at the first observation epoch are approximately:

a    = 28196.776 km  
e    = 0.7679436  
i    = 20.315 deg  
RAAN = 359.145 deg  
AOP  = 179.425 deg  
M    = 29.236 deg

These agree with the reference solution in Satellite Orbits by Montenbruck and Gill.

---

## Reusable Functions

This exercise uses the following reusable project functions:

coordinates/
    az_el_range_to_enz.py
    enz_to_ecef.py
    ecef_to_geodetic.py
    ecef_to_eci.py

time_systems/
    calendar_to_mjd.py
    mjd_to_greenwich_hour_angle.py

astrodynamics/
    gauss_sector_triangle_ratio.py
    gauss_two_position_elements.py

The exercise script itself contains only the observation data and the sequence connecting these general-purpose functions.
