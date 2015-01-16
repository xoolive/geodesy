# -*- coding: utf-8 -*-
""" Benchmarks SSE2 functions for geodesy."""

import geodesy.sphere as geo
from chrono import chronometer


@chronometer
def distance_vec(p1, p2):
    """ Do it! """
    geo.distance_fast(p1, p2)


@chronometer
def distance_nonvec():
    """ Do it! """
    [geo.distance((0., i / 60.), (0, (i + 1) / 60)) for i in range(size)]

size = 360 * 60 - 1
nbloops = 100

p1 = [(0, i / 60.) for i in range(size)]
p2 = [(0, (i + 1) / 60.) for i in range(size)]

for j in xrange(nbloops):
    distance_vec(p1, p2)
    distance_nonvec()

print "distance:\t vectorised %g regular %g speedup %f" % \
    (distance_vec.acc_tps, distance_nonvec.acc_tps,
     (distance_nonvec.acc_tps / distance_vec.acc_tps))
