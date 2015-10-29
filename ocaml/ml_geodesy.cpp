extern "C" {
#  include <caml/memory.h>
#  include <caml/alloc.h>
} // extern

#include "spherical_geodesy.h"
#include "wgs84_geodesy.h"

typedef enum _units { DEG, RAD } units;
static double pi = acos(-1.);

#define Val_none Val_int(0)
#define Some_val(v) Field(v, 0)

extern "C"
value sph_distance (
    value unit, value v_lat1, value v_lon1, value v_lat2, value v_lon2)
{
  CAMLparam5(unit, v_lat1, v_lon1, v_lat2, v_lon2);
  CAMLlocal1(v_res);

  double distance, bearing1, bearing2;
  double lat1, lon1, lat2, lon2;

  lat1 = Double_val(v_lat1); lon1 = Double_val(v_lon1);
  lat2 = Double_val(v_lat2); lon2 = Double_val(v_lon2);

  if (unit != Val_none && Int_val(Some_val(unit)) == DEG)
  {
    lat1 = lat1 * pi/180.;
    lon1 = lon1 * pi/180.;
    lat2 = lat2 * pi/180.;
    lon2 = lon2 * pi/180.;
  }

  SphericalGeodesy_d::distanceAndBearing(
      lat1, lon1, lat2, lon2,
      distance, bearing1, bearing2);

  v_res = copy_double(distance);

  CAMLreturn(v_res);
}

extern "C"
value sph_bearing (
    value unit, value v_lat1, value v_lon1, value v_lat2, value v_lon2)
{
  CAMLparam5(unit, v_lat1, v_lon1, v_lat2, v_lon2);
  CAMLlocal1(v_res);

  double distance, bearing1, bearing2;
  double lat1, lon1, lat2, lon2;

  lat1 = Double_val(v_lat1); lon1 = Double_val(v_lon1);
  lat2 = Double_val(v_lat2); lon2 = Double_val(v_lon2);

  if (unit != Val_none && Int_val(Some_val(unit)) == DEG)
  {
    lat1 = lat1 * pi/180.;
    lon1 = lon1 * pi/180.;
    lat2 = lat2 * pi/180.;
    lon2 = lon2 * pi/180.;
  }

  SphericalGeodesy_d::distanceAndBearing(
      lat1, lon1, lat2, lon2,
      distance, bearing1, bearing2);

  if (unit != Val_none && Int_val(Some_val(unit)) == DEG)
    bearing1 = bearing1 / 180. * pi;

  v_res = copy_double(bearing1);

  CAMLreturn(v_res);
}


extern "C"
value sph_distance_bearing (
    value unit, value v_lat1, value v_lon1, value v_lat2, value v_lon2)
{
  CAMLparam5(unit, v_lat1, v_lon1, v_lat2, v_lon2);
  CAMLlocal1(v_res);

  double distance, bearing1, bearing2;
  double lat1, lon1, lat2, lon2;

  lat1 = Double_val(v_lat1); lon1 = Double_val(v_lon1);
  lat2 = Double_val(v_lat2); lon2 = Double_val(v_lon2);

  if (unit != Val_none && Int_val(Some_val(unit)) == DEG)
  {
    lat1 = lat1 * pi/180.;
    lon1 = lon1 * pi/180.;
    lat2 = lat2 * pi/180.;
    lon2 = lon2 * pi/180.;
  }

  SphericalGeodesy_d::distanceAndBearing(
      lat1, lon1, lat2, lon2,
      distance, bearing1, bearing2);

  if (unit != Val_none && Int_val(Some_val(unit)) == DEG)
    bearing1 = bearing1 / 180. * pi;

  v_res = caml_alloc_tuple(2);

  Store_field(v_res, 0, copy_double(distance));
  Store_field(v_res, 1, copy_double(bearing1));

  CAMLreturn(v_res);
}

extern "C"
value sph_destination (
    value unit, value v_lat, value v_lon, value v_h, value v_d)
{
  CAMLparam5(unit, v_lat, v_lon, v_h, v_d);
  CAMLlocal1(v_res);

  double tolat, tolon, bearing;
  double lat, lon, h;

  lat = Double_val(v_lat); lon = Double_val(v_lon); h = Double_val(v_h);

  if (unit != Val_none && Int_val(Some_val(unit)) == DEG)
  {
    lat = lat * pi/180.;
    lon = lon * pi/180.;
    h   = h * pi/180.;
  }

  SphericalGeodesy_d::destination(
      lat, lon, h, Double_val(v_d),
      tolat, tolon, bearing);

  if (unit != Val_none && Int_val(Some_val(unit)) == DEG)
  {
    tolat = tolat * 180./pi;
    tolon = tolon * 180./pi;
  }

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
value wgs84_distance (
    value unit, value v_lat1, value v_lon1, value v_lat2, value v_lon2)
{
  CAMLparam5(unit, v_lat1, v_lon1, v_lat2, v_lon2);
  CAMLlocal1(v_res);

  double distance, bearing1, bearing2;
  double lat1, lon1, lat2, lon2;

  lat1 = Double_val(v_lat1); lon1 = Double_val(v_lon1);
  lat2 = Double_val(v_lat2); lon2 = Double_val(v_lon2);

  if (unit != Val_none && Int_val(Some_val(unit)) == DEG)
  {
    lat1 = lat1 * pi/180.;
    lon1 = lon1 * pi/180.;
    lat2 = lat2 * pi/180.;
    lon2 = lon2 * pi/180.;
  }

  WGS84Geodesy_d::distanceAndBearing(
      lat1, lon1, lat2, lon2,
      distance, bearing1, bearing2);

  v_res = copy_double(distance);

  CAMLreturn(v_res);
}

extern "C"
value wgs84_bearing (
    value unit, value v_lat1, value v_lon1, value v_lat2, value v_lon2)
{
  CAMLparam5(unit, v_lat1, v_lon1, v_lat2, v_lon2);
  CAMLlocal1(v_res);

  double distance, bearing1, bearing2;
  double lat1, lon1, lat2, lon2;

  lat1 = Double_val(v_lat1); lon1 = Double_val(v_lon1);
  lat2 = Double_val(v_lat2); lon2 = Double_val(v_lon2);

  if (unit != Val_none && Int_val(Some_val(unit)) == DEG)
  {
    lat1 = lat1 * pi/180.;
    lon1 = lon1 * pi/180.;
    lat2 = lat2 * pi/180.;
    lon2 = lon2 * pi/180.;
  }

  WGS84Geodesy_d::distanceAndBearing(
      lat1, lon1, lat2, lon2,
      distance, bearing1, bearing2);

  if (unit != Val_none && Int_val(Some_val(unit)) == DEG)
    bearing1 = bearing1 / 180. * pi;

  v_res = copy_double(bearing1);

  CAMLreturn(v_res);
}

extern "C"
value wgs84_distance_bearing (
    value unit, value v_lat1, value v_lon1, value v_lat2, value v_lon2)
{
  CAMLparam5(unit, v_lat1, v_lon1, v_lat2, v_lon2);
  CAMLlocal1(v_res);

  double distance, bearing1, bearing2;
  double lat1, lon1, lat2, lon2;

  lat1 = Double_val(v_lat1); lon1 = Double_val(v_lon1);
  lat2 = Double_val(v_lat2); lon2 = Double_val(v_lon2);

  if (unit != Val_none && Int_val(Some_val(unit)) == DEG)
  {
    lat1 = lat1 * pi/180.;
    lon1 = lon1 * pi/180.;
    lat2 = lat2 * pi/180.;
    lon2 = lon2 * pi/180.;
  }

  WGS84Geodesy_d::distanceAndBearing(
      lat1, lon1, lat2, lon2,
      distance, bearing1, bearing2);

  if (unit != Val_none && Int_val(Some_val(unit)) == DEG)
    bearing1 = bearing1 / 180. * pi;

  v_res = caml_alloc_tuple(2);

  Store_field(v_res, 0, copy_double(distance));
  Store_field(v_res, 1, copy_double(bearing1));

  CAMLreturn(v_res);
}

extern "C"
value wgs84_destination (
    value unit, value v_lat, value v_lon, value v_h, value v_d)
{
  CAMLparam5(unit, v_lat, v_lon, v_h, v_d);
  CAMLlocal1(v_res);

  double tolat, tolon, bearing;
  double lat, lon, h;

  lat = Double_val(v_lat); lon = Double_val(v_lon); h = Double_val(v_h);

  if (unit != Val_none && Int_val(Some_val(unit)) == DEG)
  {
    lat = lat * pi/180.;
    lon = lon * pi/180.;
    h   = h * pi/180.;
  }

  WGS84Geodesy_d::destination(
      lat, lon, h, Double_val(v_d),
      tolat, tolon, bearing);

  if (unit != Val_none && Int_val(Some_val(unit)) == DEG)
  {
    tolat = tolat * 180./pi;
    tolon = tolon * 180./pi;
  }

  v_res = caml_alloc_tuple(2);

  Store_field(v_res, 0, copy_double(tolat));
  Store_field(v_res, 1, copy_double(tolon));

  CAMLreturn(v_res);
}
