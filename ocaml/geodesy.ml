type units = | Deg | Rad

module SphericalGeodesy = struct
    external distance :
      ?units:units -> float -> float -> float -> float -> float
      = "sph_distance"
    external bearing :
      ?units:units -> float -> float -> float -> float -> float
      = "sph_bearing"
    external distance_bearing :
      ?units:units -> float -> float -> float -> float -> float * float
      = "sph_distance_bearing"
    external destination :
      ?units:units -> float -> float -> float -> float -> float * float
      = "sph_destination"
    external crosstrack :
      float -> float -> float -> float -> float -> float * float * float
      = "sph_crosstrack"
  end

module WGS84Geodesy = struct
    external distance :
      ?units:units -> float -> float -> float -> float -> float
      = "wgs84_distance"
    external bearing :
      ?units:units -> float -> float -> float -> float -> float
      = "wgs84_bearing"
    external distance_bearing :
      ?units:units -> float -> float -> float -> float -> float * float
      = "wgs84_distance_bearing"
    external destination :
      ?units:units -> float -> float -> float -> float -> float * float
      = "wgs84_destination"
  end
