extern "C" {
#  include <caml/memory.h>
#  include <caml/alloc.h>
#  include <caml/mlvalues.h>
#  include <caml/bigarray.h>
#  include <caml/fail.h>
} // extern

// This is not set by Microsoft 64-bit compiler...
#if _M_IX86_FP == 2
#  define __SSE2__
#endif

#ifdef __SSE2__
# include <xmmintrin.h>
#endif

#include "spherical_geodesy.h"
#include "wgs84_geodesy.h"

extern "C"
value sph_distance (value v_lat1, value v_lon1, value v_lat2, value v_lon2)
{
  CAMLparam4(v_lat1, v_lon1, v_lat2, v_lon2);
  CAMLlocal1(v_res);

  double distance, bearing1, bearing2;

  SphericalGeodesy_d::distanceAndBearing(
      Double_val(v_lat1), Double_val(v_lon1),
      Double_val(v_lat2), Double_val(v_lon2),
      distance, bearing1, bearing2);

  v_res = copy_double(distance);

  CAMLreturn(v_res);
}

extern "C"
value sph_bearing (value v_lat1, value v_lon1, value v_lat2, value v_lon2)
{
  CAMLparam4(v_lat1, v_lon1, v_lat2, v_lon2);
  CAMLlocal1(v_res);

  double distance, bearing1, bearing2;

  SphericalGeodesy_d::distanceAndBearing(
      Double_val(v_lat1), Double_val(v_lon1),
      Double_val(v_lat2), Double_val(v_lon2),
      distance, bearing1, bearing2);

  v_res = copy_double(bearing1);

  CAMLreturn(v_res);
}


extern "C"
value sph_distance_fast(value v_lat1, value v_lon1, value v_lat2, value v_lon2)
{
  CAMLparam4(v_lat1, v_lon1, v_lat2, v_lon2);
  CAMLlocal1(v_res);

  mlsize_t len  = Wosize_val(v_lat1) / Double_wosize;
  mlsize_t len2 = Wosize_val(v_lon1) / Double_wosize;
  mlsize_t len3 = Wosize_val(v_lat2) / Double_wosize;
  mlsize_t len4 = Wosize_val(v_lon2) / Double_wosize;

  if (len != len2 || len != len3 || len != len4)
    caml_invalid_argument("Incompatible sizes");

  size_t i = 0;
  v_res = caml_alloc(len * Double_wosize, Double_array_tag);

#ifdef __SSE2__

  __m128d distance, bearing1, bearing2;
  for (; i < len-1; i += 2)
  {
    SphericalGeodesy<__m128d>::distanceAndBearing(
        _mm_loadu_pd((double*) v_lat1 + i),
        _mm_loadu_pd((double*) v_lon1 + i),
        _mm_loadu_pd((double*) v_lat2 + i),
        _mm_loadu_pd((double*) v_lon2 + i),
        distance, bearing1, bearing2);
    _mm_storeu_pd((double*) v_res + i, distance);
  }

#endif

  double b1, b2, res;

  for ( ; i < len; ++i)
  {
    SphericalGeodesy_d::distanceAndBearing(
        Double_field(v_lat1, i), Double_field(v_lon1, i),
        Double_field(v_lat2, i), Double_field(v_lon2, i),
        res, b1, b2);
    Store_double_field(v_res, i, res);
  }

  CAMLreturn(v_res);
}

extern "C"
value sph_destination (value v_lat, value v_lon, value v_h, value v_d)
{
  CAMLparam4(v_lat, v_lon, v_h, v_d);
  CAMLlocal1(v_res);

  double tolat, tolon, bearing;
  SphericalGeodesy_d::destination(
      Double_val(v_lat), Double_val(v_lon),
      Double_val(v_h)  , Double_val(v_d),
      tolat, tolon, bearing);

  v_res = caml_alloc_tuple(2);

  Store_field(v_res, 0, copy_double(tolat));
  Store_field(v_res, 1, copy_double(tolon));

  CAMLreturn(v_res);
}

extern "C"
value sph_destination_fast(value v_lat, value v_lon, value v_h, value v_d)
{
  CAMLparam4(v_lat, v_lon, v_h, v_d);
  CAMLlocal1(v_res);

  mlsize_t len  = Wosize_val(v_lat) / Double_wosize;
  mlsize_t len2 = Wosize_val(v_lon) / Double_wosize;
  mlsize_t len3 = Wosize_val(v_h) / Double_wosize;
  mlsize_t len4 = Wosize_val(v_d) / Double_wosize;

  if (len != len2 || len != len3 || len != len4)
    caml_invalid_argument("Incompatible sizes");
  value v_reslat, v_reslon;

  v_reslat = caml_alloc(len * Double_wosize, Double_array_tag);
  v_reslon = caml_alloc(len * Double_wosize, Double_array_tag);

  size_t i = 0;

#ifdef __SSE2__

  __m128d lat2, lon2, bearing2;
  for ( ; i < len-1; i += 2)
  {
    SphericalGeodesy<__m128d>::destination(
        _mm_loadu_pd((double*) v_lat + i),
        _mm_loadu_pd((double*) v_lon + i),
        _mm_loadu_pd((double*) v_h + i),
        _mm_loadu_pd((double*) v_d + i),
        lat2, lon2, bearing2);
    _mm_storeu_pd((double*) v_reslat + i, lat2);
    _mm_storeu_pd((double*) v_reslon + i, lon2);
  }

#endif

  double res_lat, res_lon, b;

  for (; i < len; ++i)
  {
    SphericalGeodesy_d::destination(
        Double_field(v_lat, i), Double_field(v_lon, i),
        Double_field(v_h, i), Double_field(v_d, i),
        res_lat, res_lon, b);
    Store_double_field(v_reslat, i, res_lat);
    Store_double_field(v_reslon, i, res_lon);
  }


  v_res = caml_alloc_tuple(2);

  Store_field(v_res, 0, v_reslat);
  Store_field(v_res, 1, v_reslon);

  CAMLreturn(v_res);
}

extern "C"
value sph_crosstrack(value v_lat1, value v_lon1,
                     value v_latr, value v_lonr, value v_brng)
{
  CAMLparam5(v_lat1, v_lon1, v_latr, v_lonr, v_brng);
  CAMLlocal1(v_res);

  double tolat, tolon, distance;
  SphericalGeodesy_d::crosstrack(
      Double_val(v_lat1), Double_val(v_lon1),
      Double_val(v_latr), Double_val(v_lonr), Double_val(v_brng),
      tolat, tolon, distance);

  v_res = caml_alloc_tuple(3);

  Store_field(v_res, 0, copy_double(tolat));
  Store_field(v_res, 1, copy_double(tolon));
  Store_field(v_res, 2, copy_double(distance));

  CAMLreturn(v_res);

}

extern "C"
value wgs84_distance (value v_lat1, value v_lon1, value v_lat2, value v_lon2)
{
  CAMLparam4(v_lat1, v_lon1, v_lat2, v_lon2);
  CAMLlocal1(v_res);

  double distance;
  double bearing1;
  double bearing2;

  WGS84Geodesy_d::distanceAndBearing(
      Double_val(v_lat1), Double_val(v_lon1),
      Double_val(v_lat2), Double_val(v_lon2),
      distance, bearing1, bearing2);

  v_res = copy_double(distance);

  CAMLreturn(v_res);
}

extern "C"
value wgs84_destination (value v_lat, value v_lon, value v_h, value v_d)
{
  CAMLparam4(v_lat, v_lon, v_h, v_d);
  CAMLlocal1(v_res);

  double tolat;
  double tolon;
  double bearing;

  WGS84Geodesy_d::destination(
      Double_val(v_lat), Double_val(v_lon),
      Double_val(v_h), Double_val(v_d),
      tolat, tolon, bearing);

  v_res = caml_alloc_tuple(2);

  Store_field(v_res, 0, copy_double(tolat));
  Store_field(v_res, 1, copy_double(tolon));

  CAMLreturn(v_res);
}
