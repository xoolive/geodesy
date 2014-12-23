extern "C" {
#  include <caml/memory.h>
#  include <caml/alloc.h>
#  include <caml/mlvalues.h>
#  include <caml/bigarray.h>
#  include <caml/fail.h>
} // extern

#include "spherical_geodesy.h"
#include "wgs84_geodesy.h"

// This is not set by Microsoft 64-bit compiler...
#if _M_IX86_FP == 2
#  define __SSE2__
#endif

#ifdef __SSE2__
# include <xmmintrin.h>
#endif

extern "C"
value sph_distance (value v_lat1, value v_lon1, value v_lat2, value v_lon2)
{
  CAMLparam4(v_lat1, v_lon1, v_lat2, v_lon2);
  CAMLlocal1(v_res);

  float distance, bearing1, bearing2;

  SphericalGeodesy_f::distanceAndBearing(
      (float) Double_val(v_lat1), (float) Double_val(v_lon1),
      (float) Double_val(v_lat2), (float) Double_val(v_lon2),
      distance, bearing1, bearing2);

  v_res = copy_double(distance);

  CAMLreturn(v_res);
}

extern "C"
value sph_bearing (value v_lat1, value v_lon1, value v_lat2, value v_lon2)
{
  CAMLparam4(v_lat1, v_lon1, v_lat2, v_lon2);
  CAMLlocal1(v_res);

  float distance, bearing1, bearing2;

  SphericalGeodesy_f::distanceAndBearing(
      (float) Double_val(v_lat1), (float) Double_val(v_lon1),
      (float) Double_val(v_lat2), (float) Double_val(v_lon2),
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

#ifdef __SSE2__
  {
    float* lat1 = (float*) _mm_malloc(len * sizeof(float), 16);
    float* lon1 = (float*) _mm_malloc(len * sizeof(float), 16);
    float* lat2 = (float*) _mm_malloc(len * sizeof(float), 16);
    float* lon2 = (float*) _mm_malloc(len * sizeof(float), 16);
    float* res = (float*) _mm_malloc(len * sizeof(float), 16);

    for (size_t i = 0; i < len; ++i)
    {
      lat1[i] = (float) Double_field(v_lat1, i);
      lon1[i] = (float) Double_field(v_lon1, i);
      lat2[i] = (float) Double_field(v_lat2, i);
      lon2[i] = (float) Double_field(v_lon2, i);
    }

    __m128 distance, bearing1, bearing2;
    size_t i;
    for (i = 0; i < len-3; i += 4)
    {
      SphericalGeodesy<__m128>::distanceAndBearing(
          _mm_loadu_ps(lat1 + i), _mm_loadu_ps(lon1 + i),
          _mm_loadu_ps(lat2 + i), _mm_loadu_ps(lon2 + i),
          distance, bearing1, bearing2);
      _mm_storeu_ps(res + i, distance);
    }

    float b1, b2;
    for ( ; i<len; ++i)
    {
      SphericalGeodesy_f::distanceAndBearing(
          lat1[i], lon1[i], lat2[i], lon2[i], res[i], b1, b2);
    }

    v_res = caml_alloc(len * Double_wosize, Double_array_tag);

    for (size_t i = 0; i < len; ++i)
      Store_double_field(v_res, i, res[i]);

    _mm_free(lat1); _mm_free(lon1);
    _mm_free(lat2); _mm_free(lon2);
    _mm_free(res);
  }

#else

  {
    float b1, b2, res;
    v_res = caml_alloc(len * Double_wosize, Double_array_tag);

    for (size_t i = 0; i < len; ++i)
    {
      SphericalGeodesy_f::distanceAndBearing(
          (float) Double_field(v_lat1, i), (float) Double_field(v_lon1, i),
          (float) Double_field(v_lat2, i), (float) Double_field(v_lon2, i),
          res, b1, b2);
      Store_double_field(v_res, i, res);
    }
  }

#endif

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

#ifdef __SSE2__

  {
    float* lat     = (float*) _mm_malloc(len * sizeof(float), 16);
    float* lon     = (float*) _mm_malloc(len * sizeof(float), 16);
    float* h       = (float*) _mm_malloc(len * sizeof(float), 16);
    float* d       = (float*) _mm_malloc(len * sizeof(float), 16);
    float* res_lat = (float*) _mm_malloc(len * sizeof(float), 16);
    float* res_lon = (float*) _mm_malloc(len * sizeof(float), 16);

    for (size_t i = 0; i < len; ++i)
    {
      lat[i] = (float) Double_field(v_lat, i);
      lon[i] = (float) Double_field(v_lon, i);
      h[i]   = (float) Double_field(v_h, i);
      d[i]   = (float) Double_field(v_d, i);
    }

    __m128 lat2, lon2, bearing2;
    size_t i;
    for (i = 0; i < len-3; i += 4)
    {
      SphericalGeodesy<__m128>::destination(
          _mm_loadu_ps(lat + i), _mm_loadu_ps(lon + i),
          _mm_loadu_ps(h + i), _mm_loadu_ps(d + i),
          lat2, lon2, bearing2);
      _mm_storeu_ps(res_lat + i, lat2);
      _mm_storeu_ps(res_lon + i, lon2);
    }

    float b2;
    for ( ; i<len; ++i)
    {
      SphericalGeodesy_f::destination(
          lat[i], lon[i], h[i], d[i], res_lat[i], res_lon[i], b2);
    }

    for (size_t i = 0; i < len; ++i)
    {
      Store_double_field(v_reslat, i, res_lat[i]);
      Store_double_field(v_reslon, i, res_lon[i]);
    }

    _mm_free(lat); _mm_free(lon);
    _mm_free(h); _mm_free(d);
    _mm_free(res_lat); _mm_free(res_lon);
  }

#else

  {
    float res_lat, res_lon, b;

    for (size_t i = 0; i < len; ++i)
    {
      SphericalGeodesy_f::destination(
          (float) Double_field(v_lat, i), (float) Double_field(v_lon, i),
          (float) Double_field(v_h, i), (float) Double_field(v_d, i),
          res_lat, res_lon, b);
      Store_double_field(v_reslat, i, res_lat);
      Store_double_field(v_reslon, i, res_lon);
    }
  }

#endif

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
