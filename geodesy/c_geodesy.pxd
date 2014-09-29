cdef extern from "../src/spherical_geodesy.h":
    void sph_distance "SphericalGeodesy_d::distanceAndBearing"\
        (double, double, double, double, double, double, double)
    void sph_distance_f "SphericalGeodesy_f::distanceAndBearing"\
        (float, float, float, float, float, float, float)
    void sph_destination "SphericalGeodesy_d::destination"\
        (double, double, double, double, double, double, double)
    void sph_crosstrack "SphericalGeodesy_d::crosstrack"\
        (double, double, double, double,  double, double, double, double)
    void sph_intersection "SphericalGeodesy_d::intersection"\
        (double, double, double, double,  double, double, double, double)\
        except +

cdef extern from "../src/py_geodesy.h":
    void sph_distance_fast "sph_distanceAndBearing_fast"\
        (float*, float*, float*, float*, float*, float*, float*, unsigned)

cdef extern from "../src/wgs84_geodesy.h":
    void wgs84_distance "WGS84Geodesy_d::distanceAndBearing"\
        (double, double, double, double, double, double, double)
    void wgs84_destination "WGS84Geodesy_d::destination"\
        (double, double, double, double, double, double, double)

