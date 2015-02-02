#ifndef wgs84_geodesy
#define wgs84_geodesy


template <typename T>
class WGS84Geodesy
{
public:

    static void distanceAndBearing(
        const T& fromlat, const T& fromlon, const T& tolat, const T& tolon,
        T& distance, T& bearing1, T& bearing2);

    static void distanceAndBearing(
        const T* fromlat, const T* fromlon, const T* tolat, const T* tolon,
        T* distance, T* bearing1, T* bearing2, long len);

    static void destination(
        const T& fromlat, const T& fromlon, const T& bearing1, const T& dist,
        T& tolat, T& tolon, T& bearing2);

private:
    static T EarthRadius;
    static T wgs84_a;
    static T wgs84_b;
    static T wgs84_f;
    static T wgs84_w;

};

#include "wgs84_geodesy.hpp"

template<> float WGS84Geodesy<float>::EarthRadius = 6371000.0;
template<> double WGS84Geodesy<double>::EarthRadius = 6371000.0;

template<> float WGS84Geodesy<float>::wgs84_a = 6378137.;
template<> double WGS84Geodesy<double>::wgs84_a = 6378137.;

// Note that f = (a - b) / a
template<> float WGS84Geodesy<float>::wgs84_b = 6356752.3142f;
template<> double WGS84Geodesy<double>::wgs84_b = 6356752.3142;

// 1/298.257223563;
template<> float WGS84Geodesy<float>::wgs84_f = 0.00335281066474748f;
template<> double WGS84Geodesy<double>::wgs84_f = 0.00335281066474748;

// Rotation speed (7.292115e-5 rad/s)
template<> float WGS84Geodesy<float>::wgs84_w = 0.00007292115f;
template<> double WGS84Geodesy<double>::wgs84_w = 0.00007292115;

typedef WGS84Geodesy<float> WGS84Geodesy_f;
typedef WGS84Geodesy<double> WGS84Geodesy_d;

template class WGS84Geodesy<float>;
template class WGS84Geodesy<double>;

#endif

