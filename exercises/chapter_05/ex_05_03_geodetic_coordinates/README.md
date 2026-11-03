# Exercise 5.3 — Geodetic Coordinates

## Objective

Convert the Cartesian coordinates of the Diego Garcia GPS receiver
to geodetic longitude, latitude, and ellipsoidal height using the
WGS84 reference ellipsoid.

Reference: Montenbruck and Gill, Satellite Orbits,
Exercise 5.3, page 192.

## Inputs

Receiver coordinates in WGS84(G873), at coordinate epoch 1997.0:

    x =  1917032.190 m
    y =  6029782.349 m
    z =  -801376.113 m

WGS84 ellipsoid parameters from Table 5.4:

    Equatorial radius = 6378.137 km
    Flattening = 1 / 298.257223563

The input coordinates are converted from metres to kilometres
before calling the existing conversion function.

## Reusable Function

The exercise reuses:

    src/satellite_orbits/coordinates/ecef_to_geodetic.py

Interface:

    longitude, latitude, height = ecef_to_geodetic(
        r_ecef,
        earth_radius,
        flattening,
        tolerance=1e-12,
    )

Inputs:

- r_ecef: Earth-fixed Cartesian position, shape (3,), in km.
- earth_radius: Ellipsoid equatorial radius in km.
- flattening: Dimensionless ellipsoid flattening.
- tolerance: Convergence tolerance for the delta-z iteration in km.

Outputs, in order:

- Longitude in radians, positive east.
- Geodetic latitude in radians, positive north.
- Ellipsoidal height in km.

The coordinate frame and reference ellipsoid must be compatible.
For display, the solution converts the angles to degrees and
the height to metres.

No new coordinate-conversion implementation is introduced.

## Method

The existing function follows the iterative delta-z method in
Section 5.5, equations 5.85-5.88.

The ellipsoid eccentricity squared is:

    e_squared = flattening * (2 - flattening)

The iteration determines geodetic latitude and the corresponding
ellipsoid curvature radius. Longitude is obtained using atan2(y, x).

Geodetic latitude describes the direction of the ellipsoid normal;
it is not generally equal to geocentric latitude.

## Reference Results

    East longitude = +72.36312094 deg
    Latitude       =  -7.26654999 deg
    Height         = -63.667 m

The negative height places the point below the WGS84 reference
ellipsoid. It does not, by itself, determine height relative to
mean sea level.

## Verification

The solution checks agreement with the rounded textbook values:

- Longitude: absolute tolerance 1e-8 deg.
- Latitude: absolute tolerance 1e-8 deg.
- Height: absolute tolerance 0.001 m.
- Relative tolerance: zero.

It also reconstructs the Cartesian position using equations
5.83-5.84 and checks that the position error norm is below 1e-6 m.

The reconstruction check tests numerical consistency. These
tolerances do not represent the accuracy of the surveyed station
coordinates.

## Scope and Limitations

- This is a coordinate-representation conversion, not an inertial
  frame transformation or orbit propagation.
- The coordinate epoch remains 1997.0.
- No station motion or datum transformation is applied.
- Height is ellipsoidal; no geoid model is used.
- The existing inverse routine has no iteration limit or input
  validation. Invalid inputs must not be supplied.
- The Cartesian origin is undefined for this conversion.
- Longitude is undefined on the polar axis; a convention would
  be required for general-purpose polar-axis handling.


## References

Montenbruck and Gill, Satellite Orbits:

- Section 5.5: Geodetic datums.
- Equations 5.83-5.84: Geodetic-to-Cartesian reconstruction.
- Equations 5.85-5.88: Iterative Cartesian-to-geodetic conversion.
- Table 5.4, page 189: WGS84 ellipsoid parameters.
- Exercise 5.3, page 192: Station coordinates and reference result.
