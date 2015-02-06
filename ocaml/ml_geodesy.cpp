extern "C" {
#  include <caml/memory.h>
#  include <caml/alloc.h>
} // extern

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
