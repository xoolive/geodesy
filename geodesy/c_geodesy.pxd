cdef extern from "../src/spherical_geodesy.h":

    void sph_distance "SphericalGeodesy_d::distanceAndBearing"\
        (double, double, double, double, double, double, double)
    void sph_distance_vec_d "SphericalGeodesy_d::distanceAndBearing"\
        (double*, double*, double*, double*,\
         double*, double*, double*, long)
    void sph_distance_vec_f "SphericalGeodesy_f::distanceAndBearing"\
        (float*, float*, float*, float*,\
         float*, float*, float*, long)

    void sph_destination "SphericalGeodesy_d::destination"\
        (double, double, double, double, double, double, double)
    void sph_crosstrack "SphericalGeodesy_d::crosstrack"\
        (double, double, double, double,  double, double, double, double)
    void sph_intersection "SphericalGeodesy_d::intersection"\
        (double, double, double, double,  double, double, double, double)\
        except +

cdef extern from "../src/wgs84_geodesy.h":
    void wgs84_distance "WGS84Geodesy_d::distanceAndBearing"\
        (double, double, double, double, double, double, double)
    void wgs84_distance_vec_d "WGS84Geodesy_d::distanceAndBearing"\
        (double*, double*, double*, double*,\
         double*, double*, double*, long)
    void wgs84_destination "WGS84Geodesy_d::destination"\
        (double, double, double, double, double, double, double)

