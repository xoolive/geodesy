import geodesy.sphere as geo
import numpy as np
import time

size = 360 * 60 - 1

p1 = [0 for i in range(size)]
p2 = [i / 60. for i in range(size)]
p3 = [0 for i in range(size)]
p4 = [(i + 1) / 60. for i in range(size)]

p1_64 = np.array(p1, np.float64)
p2_64 = np.array(p2, np.float64)
p3_64 = np.array(p3, np.float64)
p4_64 = np.array(p4, np.float64)

p1_32 = np.array(p1, np.float32)
p2_32 = np.array(p2, np.float32)
p3_32 = np.array(p3, np.float32)
p4_32 = np.array(p4, np.float32)

t1 = time.time()

for j in range(100):
    for i in range(size):
        d = geo.distance(p1[i], p2[i], p3[i], p4[i])

t2 = time.time()

for j in range(100):
    d = geo.distance(p1, p2, p3, p4)

t3 = time.time()

for j in range(100):
    d = geo.distance(p1_64, p2_64, p3_64, p4_64)

t4 = time.time()

for j in range(100):
    d = geo.distance(p1_32, p2_32, p3_32, p4_32)

t5 = time.time()

for j in range(100):
    d = geo.distance(p1_32, p2_64, p3_32, p4_64)

t6 = time.time()


print ("iterate in Python                              %.2f  sec" % (t2-t1))
print ("iterate in geodesy                             %.2f  sec" % (t3-t2))
print ("passing numpy arrays (double precision)        %.3f sec" % (t4-t3))
print ("passing numpy arrays (simple precision)        %.3f sec" % (t5-t4))
print ("downcasting numpy arrays to simple precision   %.3f sec" % (t6-t5))
