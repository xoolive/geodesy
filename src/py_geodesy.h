#ifndef py_geodesy
#define py_geodesy

#include "spherical_geodesy.h"

void sph_distanceAndBearing_fast(
    float* lat1, float* lon1, float* lat2, float* lon2,
    float* distance, float* b1, float* b2, unsigned int length)
{

  size_t i = 0;

#ifdef __SSE2__
  __m128 d, f1, f2;
  for (; i < length; i += 4)
  {
    SphericalGeodesy<__m128>::distanceAndBearing(
        _mm_load_ps(lat1+i), _mm_load_ps(lon1+i),
        _mm_load_ps(lat2+i), _mm_load_ps(lon2+i), d, f1, f2);
    _mm_storeu_ps(distance + i, d);
  }
#endif

  for (; i < length; i += 1)
    SphericalGeodesy_f::distanceAndBearing(
        lat1[i], lon1[i], lat2[i], lon2[i], distance[i], b1[i], b2[i]);

}

#endif
