#ifndef spherical_geodesy
#define spherical_geodesy

// #  if _M_IX86_FP == 2
#  if defined(MS_WIN64)
#  define __SSE2__
#  endif

class IntersectionException {};

template <typename T>
class SphericalGeodesy
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

    static void destination(
        const T* fromlat, const T* fromlon, const T* bearing1, const T* dist,
        T* tolat, T* tolon, T* bearing2, long len);

    static void crosstrack(
        const T& fromlat, const T& fromlon,
        const T& radlat, const T& radlon, const T& radbearing,
        T& tolat, T& tolon, T& distance);

    static void intersection(
        const T& lat1, const T& lon1, const T& bearing1,
        const T& lat2, const T& lon2, const T& bearing2,
        T& lat, T& lon);

private:
    static T EarthRadius;
};

#ifdef __SSE2__
# include <emmintrin.h>
  template<> __m128 SphericalGeodesy<__m128>::EarthRadius =
    _mm_set1_ps(6371000.0);
  template<> __m128d SphericalGeodesy<__m128d>::EarthRadius =
    _mm_set1_pd(6371000.0);
#endif

#include "spherical_geodesy.hpp"

template<> float SphericalGeodesy<float>::EarthRadius = 6371000.0;
template<> double SphericalGeodesy<double>::EarthRadius = 6371000.0;

typedef SphericalGeodesy<float> SphericalGeodesy_f;
typedef SphericalGeodesy<double> SphericalGeodesy_d;

template class SphericalGeodesy<float>;
template class SphericalGeodesy<double>;

#endif

