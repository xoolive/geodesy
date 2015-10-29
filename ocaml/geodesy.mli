(** Computing with great circles for Spherical and WGS84 models. *)


(** The type for degrees and radians *)
type units = Deg | Rad

(**
   This module gives access to a series of complex formulas related to geometry
   of great circles.

   See also: [http://www.movable-type.co.uk/scripts/latlong.html]

   By default, the interface works with radians, but a [~units:Deg] may be
   positioned.
 *)
module SphericalGeodesy :
  sig
      (**
        Computes the shortest distance between two positions defined by their
        latitudes and longitudes. The method (Haversine) assumes a spherical
        Earth of radius 6,371,000 m.

        As a reference, the minute of arc in latitude at the equator defines the
         nautical mile.
       *)
    external distance :
      ?units:units -> float -> float -> float -> float -> float
      = "sph_distance"

      (**
        Computes the initial bearing (direction, azimuth) from position 1 to
        position 2 defined by their latitudes and longitudes, following a great
        circle.

        Note that unlike a loxodrome, the bearing changes continously along a
        great circle.
       *)
    external bearing :
      ?units:units -> float -> float -> float -> float -> float
      = "sph_bearing"

      (**
        Computes simultaneously distance and bearing.

        See documentation for resp. distance and bearing.
       *)
    external distance_bearing :
      ?units:units -> float -> float -> float -> float -> float * float
      = "sph_distance_bearing"

      (**
        Computes the destination point travelling along a great circle given a
        start position (lat, lon), an initial bearing and a distance.
       *)
    external destination :
      ?units:units -> float -> float -> float -> float -> float * float
      = "sph_destination"

      (**
        Computes the cross-track distance of a point p3 to a great circle,
        defined by two positions p1 and p2. The output contains the projected
        point and the distance.

        This function does not expose a similar flexible interface as the other
        ones. (I was lazy)
       *)
    external crosstrack :
      float -> float -> float -> float -> float -> float * float * float
      = "sph_crosstrack"
  end

(**
  This module computes basic geometry operations on the WGS84 model of Earth
  Ellipsoid.

  See also: [http://www.movable-type.co.uk/scripts/latlong-vincenty.html] and
  [http://www.movable-type.co.uk/scripts/latlong-vincenty-direct.html]

  By default, the interface works with radians, but a [~units:Deg] may be
  positioned.
 *)
module WGS84Geodesy :
  sig
      (**
        Computes the shortest distance between two positions defined by their
        latitudes and longitudes. The method (Vincenty) assumes a WGS84
        Ellipsoidal Earth.

        As a reference, the minute of arc in latitude at the equator defines
        the nautical mile.
       *)
    external distance :
      ?units:units -> float -> float -> float -> float -> float
      = "wgs84_distance"

      (**
        Computes the initial bearing (direction, azimuth) from position 1 to
        position 2 defined by their latitudes and longitudes, following a great
        circle. Note the bearing changes continously along a great circle. The
        method (Vincenty) assumes a WGS84 Ellipsoidal Earth.
       *)
    external bearing :
      ?units:units -> float -> float -> float -> float -> float
      = "wgs84_bearing"

      (**
        Computes simultaneously distance and bearing.

        See documentation for resp. distance and bearing.
       *)
    external distance_bearing :
      ?units:units -> float -> float -> float -> float -> float * float
      = "wgs84_distance_bearing"

      (**
        Computes the destination point given a start position (lat, lon), an
        initial bearing and a distance using Vincenty inverse formula for
        ellipsoids.
       *)
    external destination :
      ?units:units -> float -> float -> float -> float -> float * float
      = "wgs84_destination"
  end
