# -*- coding: utf-8 -*-

""" chrono docstring.  """
import time
import functools


class chronometer(object):

    """ Chronometer decorator for functions and methods.

    Usage :
    @chronometer
    def func( any parameter ) :
      Some documentations and codes...

    # Calling the function
    func( some parameters )
    print "Time spended in function : %g" % func.acc_tps
    print "Mean time spended un the function : %g" % func.meanTimePerCall()
    func.chrono_reset()

    """

    def __init__(self, f):
        self.func = f
        self.ismeth = False
        self.acc_tps = 0.
        self.nb_call = 0

    def __call__(self, *args, **kwargs):
        t1 = time.time()
        ret = self.func(*args, **kwargs)
        t2 = time.time()
        self.acc_tps += t2 - t1
        self.nb_call += 1
        return ret

    def chrono_reset(self):
        """ Reset the chronometer for a new measurement. """
        self.nb_call = 0
        self.acc_tps = 0.

    def __repr__(self):
        """ Return documentation. """
        return self.func.__doc__

    def __get__(self, obj, objtype):
        """ Return self with instances methods. """
        if not self.ismeth:
            self.func = functools.partial(self.func, obj)
            self.ismeth = True
        return self

    def meanTimePerCall(self):
        """ Return the time spent.  """
        return self.acc_tps / self.nb_call
