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

    let size = 360 * 60 - 1 in
    let zeros = Array.make size 0. in
    let lon1 = Array.init size
                 (fun i -> radians ((float_of_int i) /. 60.)) in
    let lon2 = Array.init size
                 (fun i -> radians ((float_of_int (i+1)) /. 60.)) in

    let nm = SphericalGeodesy.distance 0. 0. 0. (radians 1./.60.) in
    assert_float nm 1853.24878 "SphericalGeodesy.distance";

    let res = SphericalGeodesy.distance_fast zeros lon1 zeros lon2 in
    let d = Array.mapi (
      fun i p -> SphericalGeodesy.distance 0. p 0. lon2.(i)
    ) lon1 in

    Array.iteri (
      fun i p -> assert_float p d.(i)
                   ("SphericalGeodesy.distance_fast " ^ (string_of_int i))
    ) res;

    let test_fail_sizes () =
      let _ = SphericalGeodesy.distance_fast (Array.make (size+1) 0.)
                  zeros zeros zeros in
       () in
    assert_raises (Invalid_argument "Incompatible sizes") test_fail_sizes;

    let nm = WGS84Geodesy.distance 0. 0. 0. (radians 1./.60.) in
    assert_float nm 1855.324846 "WGS84Geodesy.distance";

  end

let test_geodesy_destination () =
  begin

    let lat, lon =
      SphericalGeodesy.destination 0. 0. (radians 90.) (60. *. 1853.25) in

    (* weird behaviour around 0. *)
    assert_float (1. +. lat) 1. "SphericalGeodesy.destination : lat";
    assert_float (degrees lon) 1. "SphericalGeodesy.destination : lon";

    let size = 360 in
    let ones = Array.make size 1. in
    let lon = Array.init size (fun i -> radians (float_of_int i)) in
    let d = Array.make size (60. *. 1853.24878) in
    let v = Array.make size (radians 90.) in

    let rlat1, rlon1 = SphericalGeodesy.destination_fast ones lon v d in

    let la, lo = SphericalGeodesy.destination 1. lon.(2) v.(2) d.(2) in
    print_newline ();

    Array.iteri (
      fun i p ->
        assert_float p
          (snd (SphericalGeodesy.destination 1. lon.(i) v.(i) d.(i)))
          ("SphericalGeodesy.destination_fast lon " ^ (string_of_int i))
    ) rlon1;

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
