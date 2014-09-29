# distutils: language = c++
# -*- coding: utf-8 -*-

"""
    This module computes basic operations on the WGS84 model of Earth Ellipsoid.
    Although the module could have been designed to accept different models'
    parameters such as Airy (1830), International (1924), better fits
    respectivaly in the UK and in much Europe, this has not been implemented to
    avoid mistakes on datums. Usually, coordinates are given for WGS84 (which
    specifies its own datum); Greenwich meridian (aka Airy transit circle) is
    about 100m far from the zero given by GPS devices (calibrated on WGS84).

    See also: http://www.movable-type.co.uk/scripts/latlong-vincenty.html
       http://www.movable-type.co.uk/scripts/latlong-vincenty-direct.html

"""

from numpy import degrees, radians

from c_geodesy cimport wgs84_distance, wgs84_destination

def distance(p1, p2):
    """ distance(p1 [deg, deg], p2 [deg, deg]) -> [m]
        Computes the shortest distance between two positions defined by
        their latitudes and longitudes. The method (Vincenty) assumes a
        WGS84 Ellipsoidal Earth.

        As a reference, the minute of arc in latitude at the equator defines the
        nautical mile.
        >>> distance((0, 0), (0, 1./60))
        1855.3248463069199
    """
    cdef double d = 0, b1 = 0, b2 = 0
    cdef double lat1 = radians(p1[0]), lon1 = radians(p1[1])
    cdef double lat2 = radians(p2[0]), lon2 = radians(p2[1])
    wgs84_distance(lat1, lon1, lat2, lon2, d, b1, b2)
    return d

def bearing(p1, p2):
    """ bearing(p1 [deg, deg], p2 [deg, deg]) -> [deg]
        Computes the initial bearing (direction, azimuth) from position 1 to
        position 2 defined by their latitudes and longitudes, following a great
        circle. Note the bearing changes continously along a great circle. The
        method (Vincenty) assumes a WGS84 Ellipsoidal Earth.

        The following example shows the direction to the East.
        >>> bearing((0, 0), (0, 1./60))
        90.0
    """
    cdef double d = 0, b1 = 0, b2 = 0
    cdef double lat1 = radians(p1[0]), lon1 = radians(p1[1])
    cdef double lat2 = radians(p2[0]), lon2 = radians(p2[1])
    wgs84_distance(lat1, lon1, lat2, lon2, d, b1, b2)
    return degrees(b1)

def destination(p, double bearing, double distance):
    """ destination(p [deg, deg], bearing [deg], distance [m])
                    -> ([deg, deg])
        Computes the destination point given a start position (lat, lon), an
        initial bearing and a distance using Vincenty inverse formula for
        ellipsoids.

        60 nm on the Equator correspond to one degree.
        >>> p = destination((0, 0), 90, 60 * 1855.32)
        >>> "%.5f %.5f" % p
        '0.00000 1.00000'
    """
    cdef double lat = radians(p[0]), lon = radians(p[1])
    cdef double tolat = 0, tolon = 0, tmp = 0
    wgs84_destination(lat, lon, radians(bearing), distance, tolat, tolon, tmp)
    return (degrees(tolat), degrees(tolon))

