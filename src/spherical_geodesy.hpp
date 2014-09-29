#include <cmath>
#include <iostream>

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

  const T sin_dlat_2 = sin(dlat / 2.0);
  const T sin_dlon_2 = sin(dlon / 2.0);
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

#endif

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
#endif

template<typename T>
void SphericalGeodesy<T>::crosstrack(
    const T& fromlat, const T& fromlon,
    const T& rlat, const T& rlon, const T& rbrng,
    T& tolat, T& tolon, T& distance)
{

  static T PI = acos(-1.);

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

  T dist12 = 2.*asin(sqrt(sin(dlat/2.) * sin(dlat/2.) +
                          cos(lat1) * cos(lat2) *
                          sin(dlon/2.) * sin(dlon/2.)));

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
  if (sin(alpha1)*sin(alpha2) < 0) throw IntersectionException();

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


