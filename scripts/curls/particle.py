



import numpy as np
import math


filename = "curl_0.particle"

file = open(filename,'r')


x = []
y = []
z = []
r = []

for line in file.readlines():
    values = line.split()
    x.append(float(values[0]))
    y.append(float(values[1]))
    z.append(float(values[2]))
    r.append(float(values[3]))


file.close()



x = np.array(x)
y = np.array(y)
z = np.array(z)
r = np.array(r)

x = x - np.mean(x)
y = y - np.mean(y)
z = z - np.mean(z)

# for (xx,yy,zz,rr) in zip(x,y,z,r):
#     print('{:.2e}'.format(xx),'{:.2e}'.format(yy),'{:.2e}'.format(zz),'{:.2e}'.format(rr))

print("AABB box lengths")
print(np.max(x)-np.min(x))
print(np.max(y)-np.min(y))
print(np.max(z)-np.min(z))

max_dis = 0
max_xcm_dis = 0


xcm = np.mean(x)
ycm = np.mean(y)
zcm = np.mean(z)

for (xx1,yy1,zz1,rr1) in zip(x,y,z,r):
    cm_distance = math.sqrt((xx1-xcm)**2 + (yy1-ycm)**2 + (zz1-zcm)**2)
    max_xcm_dis = max(cm_distance,max_xcm_dis)

    for (xx2,yy2,zz2,rr2) in zip(x,y,z,r):
        distance = math.sqrt((xx1-xx2)**2 + (yy1-yy2)**2 + (zz1-zz2)**2)
        max_dis = max(distance,max_dis)

print("max distance from any two atoms")
print(max_dis)
print("max distance from xcm to furthest atom")
print(max_xcm_dis)

print("expected min bin size")
print(max(np.max(r) * 2 * 1.1, max_xcm_dis * 1.1))


min_radius = np.min(r)


min_radius = 3e-05

poisson_ratio = 0.3
density = 2650.000
E = 3.70000e+09
alphap = 0.1631*poisson_ratio + 0.876605


G = E / (2.0 * (1.0 + poisson_ratio))
min_dt = math.pi * min_radius * math.sqrt(density / G) / alphap
print("max delta time")
print(min_dt)




