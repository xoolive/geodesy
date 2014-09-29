# distutils: language = c++
# -*- coding: utf-8 -*-

"""
    This module gives access to a series of complex formulas related to geometry
    of great circles.

    All positions are defined as a pair (lat, lon)

    See also: http://www.movable-type.co.uk/scripts/latlong.html
"""

import math
import numpy as np

cimport numpy as cnp
cimport cython
cimport cpython

from c_geodesy\
    cimport sph_distance, sph_destination, sph_intersection, sph_crosstrack,\
    sph_distance_fast

cdef double radians = math.acos(-1) / 180

def distance(p1, p2):
    """ distance(p1 [deg; deg], p2 [deg; deg]) -> [m]
        Computes the shortest distance between two positions defined by
        their latitudes and longitudes. The method (Haversine) assumes a
        spherical Earth of radius 6,371,000 m.

        As a reference, the minute of arc in latitude at the equator defines the
        nautical mile.
        >>> distance((0, 0), (0, 1./60))
        1853.2487774093122
    """
    cdef double d = 0, b1 = 0, b2 = 0
    cdef double lat1 = p1[0]*radians, lon1 = p1[1]*radians
    cdef double lat2 = p2[0]*radians, lon2 = p2[1]*radians
    sph_distance(lat1, lon1, lat2, lon2, d, b1, b2)
    return d

@cython.boundscheck(False) # I will be reasonable, I promise
@cython.wraparound (False) # no negative indices
def distance_fast(p1, p2):
    """ distance_fast(p1 [deg; deg], p2 [deg; deg]) -> [m]

        Computes the shortest distance between two positions defined by
        their latitudes and longitudes. The method (Haversine) assumes a
        spherical Earth of radius 6,371,000 m.

        p1 and p2 are any kind of iterable. The result is a numpy array.

        In order to vectorise efficiently, the operations are computed with
        single precision floating point values, hence the differences in all
        computed values.

        >>> size = 360 * 60 - 1
        >>> p1 = [(0, i / 60.) for i in range(size)]
        >>> p2 = [(0, (i + 1) / 60.) for i in range(size)]
        >>> d = distance_fast(p1, p2)
        >>> d[1:4]
        array([ 1853.2487793 ,  1853.24865723,  1853.24902344], dtype=float32)
    """
    cdef float f1 = 0, f2 = 0
    cdef long length, i = 0

    if (not cpython.PySequence_Check(p1)):
        raise TypeError("p1 must be iterable")
    if (not cpython.PySequence_Check(p2)):
        raise TypeError("p2 must be iterable")
    length = len(p1)
    if (length != len(p2)):
        raise IndexError("p1 and p2 must be of same length")

    cdef cnp.ndarray[cnp.float32_t, mode="c", ndim=1] lat1, lon1, lat2, lon2

    lat1 = np.empty((length), dtype=np.float32)
    lat2 = np.empty((length), dtype=np.float32)
    lon1 = np.empty((length), dtype=np.float32)
    lon2 = np.empty((length), dtype=np.float32)

    cdef float* pt_lat1 = <float*> cnp.PyArray_DATA(lat1)
    cdef float* pt_lon1 = <float*> cnp.PyArray_DATA(lon1)
    cdef float* pt_lat2 = <float*> cnp.PyArray_DATA(lat2)
    cdef float* pt_lon2 = <float*> cnp.PyArray_DATA(lon2)
    cdef object item1, item2
    cdef double x

    for i in range(length):
        item1 = cpython.PySequence_ITEM(p1, i)
        # PyFloat_AS_DOUBLE is much faster but without any type cast
        x = cpython.PyFloat_AsDouble(item1[0])
        pt_lat1[i] = x * radians
        x = cpython.PyFloat_AsDouble(item1[1])
        pt_lon1[i] = x * radians

        item2 = cpython.PySequence_ITEM(p2, i)
        x = cpython.PyFloat_AsDouble(item2[0])
        pt_lat2[i] = x * radians
        x = cpython.PyFloat_AsDouble(item2[1])
        pt_lon2[i] = x * radians

    cdef cnp.ndarray[cnp.float32_t, mode="c", ndim=1 ] distance
    distance = np.empty(length, np.float32)
    b1 = np.empty(length, np.float32)
    b2 = np.empty(length, np.float32)
    cdef float* pt_distance = <float*> cnp.PyArray_DATA(distance)
    cdef float* pt_b1 = <float*> cnp.PyArray_DATA(b1)
    cdef float* pt_b2 = <float*> cnp.PyArray_DATA(b2)

    sph_distance_fast(pt_lat1, pt_lon1, pt_lat2, pt_lon2,
                      pt_distance, pt_b1, pt_b2, length)

    return distance

def bearing(p1, p2):
    """ bearing(p1 [deg; deg], p2 [deg; deg]) -> [deg]
        Computes the initial bearing (direction, azimuth) from position 1 to
        position 2 defined by their latitudes and longitudes, following a great
        circle. Note that unlike a loxodrome, the bearing changes continously
        along a great circle.

        The following example shows the direction to the East.
        >>> bearing((0, 0), (0, 1./60))
        90.0
    """
    cdef double d = 0, b1 = 0, b2 = 0
    cdef double lat1 = radians * p1[0], lon1 = radians * p1[1]
    cdef double lat2 = radians * p2[0], lon2 = radians * p2[1]
    sph_distance(lat1, lon1, lat2, lon2, d, b1, b2)
    return b1 / radians

def destination(p, double bearing, float distance):
    """ destination(p [deg; deg], bearing [deg], distance [m]) -> [deg; deg]
        Computes the destination point travelling along a great circle given a
        start position (lat, lon), an initial bearing and a distance.

        60 nm on the Equator correspond to one degree.

        >>> p = destination((0, 0), 90, 60 * 1853.25)
        >>> "%.5f %.5f" % p
        '0.00000 1.00000'
    """
    cdef double lat = radians * p[0], lon = radians * p[1]
    cdef double tolat = 0, tolon = 0, tmp = 0
    sph_destination(lat, lon, radians * bearing, distance, tolat, tolon, tmp)
    return (tolat / radians, tolon / radians)

def crosstrack(p1, p2, p3):
    """ crosstrack(p1 [deg; deg], p2 [deg; deg], p3 [deg; deg])
            -> ([deg; deg], [m])
        Computes the cross-track distance of a point p3 to a great circle,
        defined by two positions p1 and p2. The output contains the projected
        point and the distance.

        >>> p = crosstrack((0,0), (0,1), (1,0))
        >>> "%.4f %.4f %.4f" % (p[0][0], p[0][1], p[1])
        '0.0000 0.0000 111194.9266'
        """
    cdef double lat1 = radians * p1[0], lon1 = radians * p1[1]
    cdef double lat3 = radians * p3[0], lon3 = radians * p3[1]
    cdef double lat4 = 0, lon4 = 0
    cdef double b = radians * bearing(p1, p2)
    cdef double d = 0
    sph_crosstrack(lat3, lon3, lat1, lon1, b, lat4, lon4, d)
    return ((lat4 / radians, lon4 / radians), d)

def intersection(p1, b1, p2, b2):
    """ intersection(p1 [deg; deg], b1 [deg], p2 [deg; deg], b2 [deg])
                     -> [deg; deg]
        Computes the intersection of two great circles, each defined by one
        position, resp. p1, p2, and one bearing, resp. b1, b2.
        Returns None if not computable.

        >>> p = intersection((0,-2), 90, (-1,0), 0)
        >>> "%.5f %.5f" % p
        '0.00000 0.00015'
    """
    cdef double lat1 = radians * p1[0], lon1 = radians * p1[1]
    cdef double lat2 = radians * p2[0], lon2 = radians * p2[1]
    cdef double br1  = radians * b1, br2 = radians * b2
    cdef double lat = 0, lon = 0
    try:
        sph_intersection(lat1, lon1, br1, lat2, lon2, br2, lat, lon)
        return (lat / radians, lon / radians)
    except Exception:
        return None

def latitudemax(p,b):
    """ latitudemax(p [deg; deg], b [deg]) -> [deg]
        Returns the maximum latitude of a great circle path, given a position
        and a bearing on the great circle. Also known as Clairaut's equation.

        >>> latitudemax((0, 0), 0)
        90.0
    """
    return (math.acos(abs(math.sin(radians * b) *
                          math.cos(radians * p[0])))) / radians

