
module SphericalGeodesy = struct
    external distance : float -> float -> float -> float -> float
      = "sph_distance"
    external bearing : float -> float -> float -> float -> float
      = "sph_bearing"
    external destination : float -> float -> float -> float -> float * float
      = "sph_destination"
    external crosstrack :
      float -> float -> float -> float -> float -> float * float * float
      = "sph_crosstrack"
  end

module WGS84Geodesy = struct
    external distance : float -> float -> float -> float -> float
      = "wgs84_distance"
    external destination : float -> float -> float -> float -> float * float
      = "wgs84_destination"
  end
