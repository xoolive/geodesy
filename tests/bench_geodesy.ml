open Geodesy

let pi = acos(-1.)
let radians x = x *. pi /. 180.
let degrees x = x *. 180. /. pi

let nbloops = 100
let size = 360 * 60 - 1

let zeros = Array.make size 0.
let lon1 = Array.init size (fun i -> radians (float_of_int i /. 60.))
let lon2 = Array.init size (fun i -> radians (float_of_int (i+1) /. 60.))

let _ =

  let t1 = Sys.time () in

  for j = 0 to nbloops do
    let _ = SphericalGeodesy.distance_fast zeros lon1 zeros lon2 in
    ()
  done;

  let t2 = Sys.time () in

  for j = 0 to nbloops do
    let _ = Array.mapi
              (fun i p -> SphericalGeodesy.distance
                            zeros.(i) lon1.(i) zeros.(i) lon2.(i)) lon1 in
    ()
  done;

  let t3 = Sys.time () in

  Printf.printf
    "distance:\tvectorised %g regular %g speedup %f\n"
    (t2 -. t1) (t3 -. t2) ((t3 -. t2) /. (t2 -. t1));
  ()
