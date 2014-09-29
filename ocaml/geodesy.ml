
module SphericalGeodesy = struct
    external distance : float -> float -> float -> float -> float
      = "sph_distance"
    external destination : float -> float -> float -> float -> float * float
      = "sph_destination"
    external crosstrack :
      float -> float -> float -> float -> float -> float * float * float
      = "sph_crosstrack"
    external distance_fast :
      float array -> float array -> float array -> float array -> float array
      = "sph_distance_fast"
    external destination_fast :
      float array ->
      float array -> float array -> float array -> float array * float array
      = "sph_destination_fast"
  end

module WGS84Geodesy = struct
    external distance : float -> float -> float -> float -> float
      = "wgs84_distance"
    external destination : float -> float -> float -> float -> float * float
      = "wgs84_destination"
  end
