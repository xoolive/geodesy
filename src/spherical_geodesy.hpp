#include <cmath>
#include <iostream>

// This is not set by Microsoft 64-bit compiler...
#if _M_IX86_FP == 2
#  define __SSE2__
#endif

#ifdef __SSE2__
#  include "sse2_math.h"
#endif

// Microsoft compiler
#  if defined(_MSC_VER)
#    define ALIGN16_BEG __declspec(align(16))
#    define ALIGN16_END
#  elif defined(__GNUC__)
#    define ALIGN16_BEG
#    define ALIGN16_END __attribute__((aligned (16)))
#  endif

#define TWO ((T) 2.0)

template<typename T>
void SphericalGeodesy<T>::distanceAndBearing(
    const T& fromlat, const T& fromlon,
    const T& tolat, const T& tolon, T& distance,
    T& bearing1, T& bearing2)
{
  const T dlat     = tolat - fromlat;
  const T dlon     = tolon - fromlon;
  const T cos_lat1 = cos(fromlat);
  const T cos_lat2 = cos(tolat);

  const T sin_dlat_2 = sin(dlat / TWO);
  const T sin_dlon_2 = sin(dlon / TWO);
  const T a = sin_dlat_2 * sin_dlat_2 +
    sin_dlon_2 * sin_dlon_2 * cos_lat1 * cos_lat2;
  const T c = 2 * atan2(sqrt(a), sqrt(1-a));
  distance  = c * EarthRadius;

  const T y1 = sin(dlon) * cos_lat2;
  const T x1 = cos_lat1 * sin(tolat) - sin(fromlat) * cos_lat2 * cos(dlon);
  bearing1   = atan2(y1, x1);

  const T y2 = - sin(dlon) * cos_lat1;
  const T x2 = cos_lat2 * sin(fromlat) - sin(tolat) * cos_lat1 * cos(dlon);
  bearing2   = atan2(y2, x2) ;
}

#ifdef __SSE2__
template<>
void SphericalGeodesy<__m128>::distanceAndBearing(
    const __m128& fromlat, const __m128& fromlon,
    const __m128& tolat, const __m128& tolon, __m128& distance,
    __m128& bearing1, __m128& bearing2)
{
  static const __m128 one_ps = _mm_set1_ps(1.0f);
  static const __m128 two_ps = _mm_set1_ps(2.0f);
  static const __m128 two_earth = _mm_mul_ps(two_ps, EarthRadius);

//   const T dlat = tolat - fromlat;
//   const T dlon = tolon - fromlon;
  __m128 dlat = _mm_sub_ps(tolat, fromlat);
  __m128 dlon = _mm_sub_ps(tolon, fromlon);

//   const T cos_lat1 = cos(fromlat);
//   const T cos_lat2 = cos(tolat);
  ALIGN16_BEG static __m128 cos_lat1 ALIGN16_END;
  ALIGN16_BEG static __m128 cos_lat2 ALIGN16_END;
  __m128 sin_lat1 = _mm_sincos_ps(&cos_lat1, fromlat);
  __m128 sin_lat2 = _mm_sincos_ps(&cos_lat2, fromlat);

  ALIGN16_BEG __m128 cos_dlon ALIGN16_END;
  __m128 sin_dlon = _mm_sincos_ps(&cos_dlon, dlon);

//   const T sin_dlat_2 = sin(dlat / 2.0);
//   const T sin_dlon_2 = sin(dlon / 2.0);
  __m128 sin_dlat_2 = _mm_sin_ps ( _mm_div_ps( dlat, two_ps ) );
  __m128 sin_dlon_2 = _mm_sin_ps ( _mm_div_ps( dlon, two_ps ) );

//   const T a = sin_dlat_2 * sin_dlat_2 +
//     sin_dlon_2 * sin_dlon_2 * cos_lat1 * cos_lat2;
//   const T c = 2 * atan2(sqrt(a), sqrt(1-a));
//   distance = c * EarthRadius;
  __m128 a = _mm_add_ps (_mm_mul_ps ( sin_dlat_2, sin_dlat_2 ),
                          _mm_mul_ps ( _mm_mul_ps ( sin_dlon_2, sin_dlon_2 ),
                                       _mm_mul_ps ( cos_lat1, cos_lat2 )));

  __m128 c = _mm_atan2_ps (_mm_sqrt_ps(a),
                           _mm_sqrt_ps ( _mm_sub_ps ( one_ps, a) ));

  distance = _mm_mul_ps (two_earth, c);

//   const T y1 = sin(dlon) * cos_lat2;
//   const T x1 =
//     cos_lat1 * sin(tolat) - sin(fromlat) * cos_lat2 * cos(dlon);
//   bearing1 = atan2(y1, x1);
  __m128 y1   = _mm_mul_ps(sin_dlon, cos_lat2);
  __m128 x1   = _mm_mul_ps(cos_lat1, sin_lat2);
  x1     = _mm_sub_ps(x1, _mm_mul_ps(_mm_mul_ps(sin_lat1, cos_lat2), cos_dlon));
  bearing1    = _mm_atan2_ps(y1, x1);

//   const T y2 = - sin(dlon) * cos_lat1;
//   const T x2 =
//     cos_lat2 * sin(fromlat) - sin(tolat) * cos_lat1 * cos(dlon);
//   bearing2 =  atan2(y2, x2) ;

  static const __m128 mone = _mm_set1_ps(-1.0f);

  __m128 y2 = _mm_mul_ps(mone, _mm_mul_ps(sin_dlon, cos_lat1));
  __m128 x2 = _mm_mul_ps(cos_lat2, sin_lat1);
  x2     = _mm_sub_ps(x2, _mm_mul_ps(_mm_mul_ps(sin_lat2, cos_lat1), cos_dlon));
  bearing2  = _mm_atan2_ps(y2, x2);

}

template<>
void SphericalGeodesy<__m128d>::distanceAndBearing(
    const __m128d& fromlat, const __m128d& fromlon,
    const __m128d& tolat, const __m128d& tolon, __m128d& distance,
    __m128d& bearing1, __m128d& bearing2)
{
  static const __m128d one_pd = _mm_set1_pd(1.0);
  static const __m128d two_pd = _mm_set1_pd(2.0);
  static const __m128d two_earth = _mm_mul_pd(two_pd, EarthRadius);

//   const T dlat = tolat - fromlat;
//   const T dlon = tolon - fromlon;
  __m128d dlat = _mm_sub_pd(tolat, fromlat);
  __m128d dlon = _mm_sub_pd(tolon, fromlon);

//   const T cos_lat1 = cos(fromlat);
//   const T cos_lat2 = cos(tolat);
  ALIGN16_BEG static __m128d cos_lat1 ALIGN16_END;
  ALIGN16_BEG static __m128d cos_lat2 ALIGN16_END;
  __m128d sin_lat1 = _mm_sincos_pd(&cos_lat1, fromlat);
  __m128d sin_lat2 = _mm_sincos_pd(&cos_lat2, fromlat);

  ALIGN16_BEG __m128d cos_dlon ALIGN16_END;
  __m128d sin_dlon = _mm_sincos_pd(&cos_dlon, dlon);

//   const T sin_dlat_2 = sin(dlat / 2.0);
//   const T sin_dlon_2 = sin(dlon / 2.0);
  __m128d sin_dlat_2 = _mm_sin_pd ( _mm_div_pd( dlat, two_pd ) );
  __m128d sin_dlon_2 = _mm_sin_pd ( _mm_div_pd( dlon, two_pd ) );

//   const T a = sin_dlat_2 * sin_dlat_2 +
//     sin_dlon_2 * sin_dlon_2 * cos_lat1 * cos_lat2;
//   const T c = 2 * atan2(sqrt(a), sqrt(1-a));
//   distance = c * EarthRadius;
  __m128d a = _mm_add_pd (_mm_mul_pd ( sin_dlat_2, sin_dlat_2 ),
                          _mm_mul_pd ( _mm_mul_pd ( sin_dlon_2, sin_dlon_2 ),
                                       _mm_mul_pd ( cos_lat1, cos_lat2 )));

  __m128d c = _mm_atan2_pd (_mm_sqrt_pd(a),
                           _mm_sqrt_pd ( _mm_sub_pd ( one_pd, a) ));

  distance = _mm_mul_pd (two_earth, c);

//   const T y1 = sin(dlon) * cos_lat2;
//   const T x1 =
//     cos_lat1 * sin(tolat) - sin(fromlat) * cos_lat2 * cos(dlon);
//   bearing1 = atan2(y1, x1);
  __m128d y1   = _mm_mul_pd(sin_dlon, cos_lat2);
  __m128d x1   = _mm_mul_pd(cos_lat1, sin_lat2);
  x1     = _mm_sub_pd(x1, _mm_mul_pd(_mm_mul_pd(sin_lat1, cos_lat2), cos_dlon));
  bearing1    = _mm_atan2_pd(y1, x1);

//   const T y2 = - sin(dlon) * cos_lat1;
//   const T x2 =
//     cos_lat2 * sin(fromlat) - sin(tolat) * cos_lat1 * cos(dlon);
//   bearing2 =  atan2(y2, x2) ;

  static const __m128d mone = _mm_set1_pd(-1.0);

  __m128d y2 = _mm_mul_pd(mone, _mm_mul_pd(sin_dlon, cos_lat1));
  __m128d x2 = _mm_mul_pd(cos_lat2, sin_lat1);
  x2     = _mm_sub_pd(x2, _mm_mul_pd(_mm_mul_pd(sin_lat2, cos_lat1), cos_dlon));
  bearing2  = _mm_atan2_pd(y2, x2);

}

#endif

template<typename T>
void SphericalGeodesy<T>::distanceAndBearing(
    const T* fromlat, const T* fromlon, const T* tolat, const T* tolon,
    T* distance, T* bearing1, T* bearing2, long len)
{
    long i = 0;
    for (; i < len ; ++i)
    {
        distanceAndBearing(fromlat[i], fromlon[i], tolat[i], tolon[i],
                           distance[i], bearing1[i], bearing2[i]);
    }
}

/*
 * DON'T DO THAT !
 * Software trigonometry is really shitty...
 *
 * template<>
 * void SphericalGeodesy<double>::distanceAndBearing(
 *     const double* fromlat, const double* fromlon,
 *     const double* tolat, const double* tolon,
 *     double* distance, double* bearing1, double* bearing2, long len)
 * {
 *     long i = 0;
 * #ifdef __SSE2__
 *     __m128d d, b1, b2;
 *     for (; i < len - 1; i += 2)
 *     {
 *         SphericalGeodesy<__m128d>::distanceAndBearing(
 *             _mm_load_pd(fromlat+i), _mm_load_pd(fromlon+i),
 *             _mm_load_pd(tolat+i), _mm_load_pd(tolon+i), d, b1, b2);
 *         _mm_storeu_pd(distance + i, d);
 *         _mm_storeu_pd(bearing1 + i, b1);
 *         _mm_storeu_pd(bearing2 + i, b2);
 *     }
 * #endif
 *     for ( ; i < len ; ++i)
 *     {
 *         distanceAndBearing(fromlat[i], fromlon[i], tolat[i], tolon[i],
 *                            distance[i], bearing1[i], bearing2[i]);
 *     }
 * }
 *
 */

template<>
void SphericalGeodesy<float>::distanceAndBearing(
    const float* fromlat, const float* fromlon,
    const float* tolat, const float* tolon,
    float* distance, float* bearing1, float* bearing2, long len)
{
    long i = 0;
#ifdef __SSE2__
    __m128 d, b1, b2;
    for (; i < len - 3; i += 4)
    {
        SphericalGeodesy<__m128>::distanceAndBearing(
            _mm_load_ps(fromlat+i), _mm_load_ps(fromlon+i),
            _mm_load_ps(tolat+i), _mm_load_ps(tolon+i), d, b1, b2);
        _mm_storeu_ps(distance + i, d);
        _mm_storeu_ps(bearing1 + i, b1);
        _mm_storeu_ps(bearing2 + i, b2);
    }
#endif
    for ( ; i < len ; ++i)
    {
        distanceAndBearing(fromlat[i], fromlon[i], tolat[i], tolon[i],
                           distance[i], bearing1[i], bearing2[i]);
    }
}

template<typename T>
void SphericalGeodesy<T>::destination(
    const T& fromlat, const T& fromlon, const T& bearing1, const T& distance,
    T& tolat, T& tolon, T& bearing2)
{
  const T sin_lat1 = sin(fromlat);
  const T cos_lat1 = cos(fromlat);
  const T d_R = distance / EarthRadius;
  const T sin_d_R = sin(d_R);
  const T cos_d_R = cos(d_R);

  tolat = asin(sin_lat1 * cos_d_R +
               cos_lat1 * sin_d_R * cos(bearing1));
  tolon = fromlon + atan2(sin(bearing1) * sin_d_R * cos_lat1,
                          cos_d_R - sin_lat1 * sin(tolat));

  const T y2 = -sin(tolat - fromlat) * cos_lat1;
  const T x2 = cos(tolat) * sin_lat1 -
    sin(tolat) * cos_lat1 * cos(tolon - fromlon);
  bearing2 = atan2(y2, x2);
}

#ifdef __SSE2__

template<>
void SphericalGeodesy<__m128>::destination(
    const __m128& fromlat, const __m128& fromlon,
    const __m128& bearing1, const __m128& distance,
    __m128& tolat, __m128& tolon, __m128& bearing2)
{
//   const T sin_lat1 = sin(fromlat);
//   const T cos_lat1 = cos(fromlat);
  ALIGN16_BEG static __m128 cos_lat1 ALIGN16_END;
  __m128 sin_lat1 = _mm_sincos_ps(&cos_lat1, fromlat);

//   const T d_R = distance / EarthRadius;
//   const T sin_d_R = sin(d_R);
//   const T cos_d_R = cos(d_R);
  ALIGN16_BEG static __m128 cos_d_R ALIGN16_END;
  __m128 d_R = _mm_div_ps(distance, EarthRadius);
  __m128 sin_d_R = _mm_sincos_ps(&cos_d_R, d_R);


  ALIGN16_BEG static __m128 cos_b1 ALIGN16_END;
  __m128 sin_b1 = _mm_sincos_ps(&cos_b1, bearing1);

//   tolat = asin(sin_lat1 * cos_d_R + cos_lat1 * sin_d_R * cos(bearing1));
  tolat = _mm_asin_ps(
      _mm_add_ps(_mm_mul_ps(sin_lat1, cos_d_R),
                 _mm_mul_ps(_mm_mul_ps(cos_lat1, sin_d_R), cos_b1)));

//   tolon = fromlon + atan2(sin(bearing1) * sin_d_R * cos_lat1,
//                           cos_d_R - sin_lat1 * sin(tolat));
  ALIGN16_BEG static __m128 cos_lat2 ALIGN16_END;
  __m128 sin_lat2 = _mm_sincos_ps(&cos_lat2, tolat);

  tolon = _mm_atan2_ps(_mm_mul_ps(sin_b1, _mm_mul_ps(sin_d_R, cos_lat1)),
                       _mm_sub_ps(cos_d_R, _mm_mul_ps(sin_lat1, sin_lat2)));
  tolon = _mm_add_ps(tolon, fromlon);

//   const T y2 = -sin(tolat - fromlat) * cos_lat1;
//   const T x2 = cos(tolat) * sin_lat1 -
//     sin(tolat) * cos_lat1 * cos(tolon - fromlon);
//   bearing2 = atan2(y2, x2);

  __m128 y2 = _mm_mul_ps(_mm_sin_ps(_mm_sub_ps(fromlat, tolat)), cos_lat1);
  __m128 x2 = _mm_sub_ps(_mm_mul_ps(cos_lat2, sin_lat1),
                         _mm_mul_ps(_mm_mul_ps(sin_lat2, cos_lat1),
                                    _mm_cos_ps(_mm_sub_ps(tolon, fromlon))));

  bearing2 = _mm_atan2_ps(y2, x2);
}

template<>
void SphericalGeodesy<__m128d>::destination(
    const __m128d& fromlat, const __m128d& fromlon,
    const __m128d& bearing1, const __m128d& distance,
    __m128d& tolat, __m128d& tolon, __m128d& bearing2)
{
//   const T sin_lat1 = sin(fromlat);
//   const T cos_lat1 = cos(fromlat);
  ALIGN16_BEG static __m128d cos_lat1 ALIGN16_END;
  __m128d sin_lat1 = _mm_sincos_pd(&cos_lat1, fromlat);

//   const T d_R = distance / EarthRadius;
//   const T sin_d_R = sin(d_R);
//   const T cos_d_R = cos(d_R);
  ALIGN16_BEG static __m128d cos_d_R ALIGN16_END;
  __m128d d_R = _mm_div_pd(distance, EarthRadius);
  __m128d sin_d_R = _mm_sincos_pd(&cos_d_R, d_R);


  ALIGN16_BEG static __m128d cos_b1 ALIGN16_END;
  __m128d sin_b1 = _mm_sincos_pd(&cos_b1, bearing1);

//   tolat = asin(sin_lat1 * cos_d_R + cos_lat1 * sin_d_R * cos(bearing1));
  tolat = _mm_asin_pd(
      _mm_add_pd(_mm_mul_pd(sin_lat1, cos_d_R),
                 _mm_mul_pd(_mm_mul_pd(cos_lat1, sin_d_R), cos_b1)));

//   tolon = fromlon + atan2(sin(bearing1) * sin_d_R * cos_lat1,
//                           cos_d_R - sin_lat1 * sin(tolat));
  ALIGN16_BEG static __m128d cos_lat2 ALIGN16_END;
  __m128d sin_lat2 = _mm_sincos_pd(&cos_lat2, tolat);

  tolon = _mm_atan2_pd(_mm_mul_pd(sin_b1, _mm_mul_pd(sin_d_R, cos_lat1)),
                       _mm_sub_pd(cos_d_R, _mm_mul_pd(sin_lat1, sin_lat2)));
  tolon = _mm_add_pd(tolon, fromlon);

//   const T y2 = -sin(tolat - fromlat) * cos_lat1;
//   const T x2 = cos(tolat) * sin_lat1 -
//     sin(tolat) * cos_lat1 * cos(tolon - fromlon);
//   bearing2 = atan2(y2, x2);

  __m128d y2 = _mm_mul_pd(_mm_sin_pd(_mm_sub_pd(fromlat, tolat)), cos_lat1);
  __m128d x2 = _mm_sub_pd(_mm_mul_pd(cos_lat2, sin_lat1),
                         _mm_mul_pd(_mm_mul_pd(sin_lat2, cos_lat1),
                                    _mm_cos_pd(_mm_sub_pd(tolon, fromlon))));

  bearing2 = _mm_atan2_pd(y2, x2);
}
#endif

template<typename T>
void SphericalGeodesy<T>::destination(
    const T* fromlat, const T* fromlon, const T* bearing1, const T* distance,
    T* tolat, T* tolon, T* bearing2, long len)
{
    long i = 0;
    for (; i < len ; ++i)
    {
        destination(fromlat[i], fromlon[i], bearing1[i], distance[i],
                    tolat[i], tolon[i], bearing2[i]);
    }
}

template<>
void SphericalGeodesy<float>::destination(
    const float* fromlat, const float* fromlon,
    const float* bearing1, const float* distance,
    float* tolat, float* tolon, float* bearing2, long len)
{
    long i = 0;
#ifdef __SSE2__
    __m128 lat, lon, b2;
    for (; i < len - 3; i += 4)
    {
        SphericalGeodesy<__m128>::destination(
            _mm_load_ps(fromlat+i), _mm_load_ps(fromlon+i),
            _mm_load_ps(bearing1+i), _mm_load_ps(distance+i),
            lat, lon, b2);
        _mm_storeu_ps(tolat + i, lat);
        _mm_storeu_ps(tolon + i, lon);
        _mm_storeu_ps(bearing2 + i, b2);
    }
#endif
    for (; i < len ; ++i)
    {
        destination(fromlat[i], fromlon[i], bearing1[i], distance[i],
                    tolat[i], tolon[i], bearing2[i]);
    }
}

template<typename T>
void SphericalGeodesy<T>::crosstrack(
    const T& fromlat, const T& fromlon,
    const T& rlat, const T& rlon, const T& rbrng,
    T& tolat, T& tolon, T& distance)
{

  static T PI = (T) acos(-1.);

  // Compute the distance to the closest point ("cross-track error")
  T distance2, tmp, bearing2;
  distanceAndBearing(rlat, rlon, fromlat, fromlon, distance2, bearing2, tmp);

  distance = EarthRadius * asin(sin(distance2 / EarthRadius) *
                                sin(bearing2 - rbrng));
  const T distance1 =
    EarthRadius * acos(cos(distance2 / EarthRadius) /
                       cos(distance / EarthRadius));

  // Compute the closest point
  T lata, lona, latb, lonb, da, db;

  destination(rlat, rlon, distance1, rbrng, lata, lona, tmp);
  destination(rlat, rlon, distance1, rbrng + PI, latb, lonb, tmp);

  distanceAndBearing(lata, lona, fromlat, fromlon, da, tmp, tmp);
  distanceAndBearing(latb, lonb, fromlat, fromlon, db, tmp, tmp);

  tolat = (da < db) ? lata : latb;
  tolon = (da < db) ? lona : lonb;

  distance = fabs(distance);
}

template<typename T>
void SphericalGeodesy<T>::intersection(
    const T& lat1, const T& lon1, const T& brng1,
    const T& lat2, const T& lon2, const T& brng2,
    T& lat3, T& lon3)
{

  const T dlat = lat2 - lat1, dlon = lon2 - lon1;

  T dist12 = TWO * asin(sqrt(sin(dlat/TWO) * sin(dlat/TWO) +
                             cos(lat1) * cos(lat2) *
                             sin(dlon/TWO) * sin(dlon/TWO)));

  if (0 == dist12) throw IntersectionException();

  const T brngA = acos((sin(lat2) - sin(lat1)*cos(dist12)) /
                       sin(dist12) * cos(lat1));
  const T brngB = acos((sin(lat1) - sin(lat2)*cos(dist12)) /
                       sin(dist12) * cos(lat2));

  T brng12, brng21;

  if (sin(dlon) > 0 ) { brng12 = brngA; brng21 = -brngB; }
  else { brng12 = -brngA; brng21 = brngB; }

  const T alpha1 = (brng1 - brng12);
  const T alpha2 = (brng21 - brng2);

  if (sin(alpha1) == 0 && sin(alpha2) == 0) throw IntersectionException();
  if (sin(alpha1) * sin(alpha2) < 0) throw IntersectionException();

  const T alpha3 = acos(-cos(alpha1) * cos(alpha2) +
                        sin(alpha1) * sin(alpha2) * cos(dist12));

  const T dist13 = atan2(sin(dist12) * sin(alpha1) * sin(alpha2),
                         cos(alpha2) + cos(alpha1) * cos(alpha3));

  lat3 = asin(sin(lat1) * cos(dist13) +
              cos(lat1) * sin(dist13) * cos(brng1));

  const T dLon13 = atan2(sin(brng1) * sin(dist13) * cos(lat1),
                         cos(dist13) - sin(lat1) * sin(lat3));

  lon3 = lon1 + dLon13;

}


