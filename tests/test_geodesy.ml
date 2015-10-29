open OUnit
open Geodesy

let pi = acos(-1.)
let radians x = x *. pi /. 180.
let degrees x = x *. 180. /. pi

let assert_float a b mymsg = assert_equal ~cmp:(cmp_float ~epsilon:1e-5)
                               ~msg:(mymsg ^
                                     " expected: " ^ string_of_float b ^
                                     " computed: " ^ string_of_float a) a b

let test_geodesy_distance () =
  begin

    let nm = SphericalGeodesy.distance 0. 0. 0. (radians 1./.60.) in
    assert_float nm 1853.248777409312 "SphericalGeodesy.distance";

    let nm = SphericalGeodesy.distance ~units:Deg 0. 0. 0. 1./.60. in
    assert_float nm 1853.248777409312 "SphericalGeodesy.distance";

    let nm = WGS84Geodesy.distance 0. 0. 0. (radians 1./.60.) in
    assert_float nm 1855.324846 "WGS84Geodesy.distance";

  end

let test_geodesy_bearing () =
  begin

    let b = SphericalGeodesy.bearing 0. 0. 0. (radians 1./.60.) in
    assert_float b pi "SphericalGeodesy.distance";

    let b = WGS84Geodesy.bearing 0. 0. 0. (radians 1./.60.) in
    assert_float b pi "WGS84Geodesy.distance";

  end

let test_geodesy_destination () =
  begin

    let lat, lon =
      SphericalGeodesy.destination 0. 0. (radians 90.) (60. *. 1853.25) in

    (* weird behaviour around 0. *)
    assert_float (1. +. lat) 1. "SphericalGeodesy.destination : lat";
    assert_float (degrees lon) 1. "SphericalGeodesy.destination : lon";

    let lat, lon =
      WGS84Geodesy.destination 0. 0. (radians 90.) (60. *. 1855.32) in

    (* weird behaviour around 0. *)
    assert_float (1. +. lat) 1. "WGS84Geodesy.destination : lat";
    assert_float (degrees lon) 1. "WGS84Geodesy.destination : lon";

  end

let suite = "Geodesy tests" >:::
            ["distance"
               >:: test_geodesy_distance;
             "destination"
               >:: test_geodesy_destination;
            ]

exception TestFailed

let _ =
  let result = run_test_tt_main suite in

  (* ugly copy paste from oUnitUtils.ml *)
  let rec was_successful =
    function
    | [] -> ()
    | RSuccess _::t
    | RSkip _::t ->
        was_successful t

    | RFailure _::_
    | RError _::_
    | RTodo _::_ ->
        raise TestFailed in

  was_successful result
