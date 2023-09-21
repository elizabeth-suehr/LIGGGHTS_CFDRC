import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt
import math

###  Superquadric code used by lebc.py
#  This may be commented out in lebc.py because of the use of scipy.special import not on HPC



def beta(x,y):
    return gamma(x) * gamma(y) / (gamma(x+y))


def product_all(v):
    return v[0]*v[1]*v[2]*v[3]*v[4]


def volume(v):
    return 2*product_all(v)*beta(v[3]/2.0 + 1, v[3])*beta(v[4]/2.0,v[4]/2.0)





def moment_of_inertia(v):
    I_xx = 0.5 * product_all(v) * ( v[1]**2 * beta(3.0/2.0*v[4],0.5*v[4])*beta(0.5*v[3],2.0*v[3]+1) + 4.0*v[2]**2*beta(0.5*v[4],0.5*v[4] + 1) * beta(3.0/2.0*v[3], v[3] + 1) )
    I_yy = 0.5 * product_all(v) * ( v[0]**2 * beta(3.0/2.0*v[4],0.5*v[4])*beta(0.5*v[3],2.0*v[3]+1) + 4.0*v[2]**2*beta(0.5*v[4],0.5*v[4] + 1) * beta(3.0/2.0*v[3], v[3] + 1) )
    I_zz = 0.5 * product_all(v) * (v[0]**2 + v[1]**2) *beta(3.0/2.0*v[4],0.5*v[4])*beta(0.5*v[3],2.0*v[3]+1)
    return [I_xx,I_yy,I_zz]



def fexp(x,p):
    return (np.sign(x)*(np.abs(x)**p))

def tens_fld(A,B,C,P,Q):
    phi, theta = np.mgrid[0:np.pi:200j, 0:2*np.pi:200j]
    x       =A*(fexp(np.sin(phi),P)) *(fexp(np.cos(theta),Q))
    y       =B*(fexp(np.sin(phi),P)) *(fexp(np.sin(theta),Q))
    z       =C*(fexp(np.cos(phi),P))
    return x , y , z 



def plot_and_save(v):
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    x,y,z = tens_fld(v[0],v[1],v[2],v[3],v[4])
    # Plot the surface.

    max_axis = max(max(v[0],v[1]),v[2])
    ax.set_xlim(-max_axis*1.1, max_axis*1.1)
    ax.set_ylim(-max_axis*1.1, max_axis*1.1)
    ax.set_zlim(-max_axis*1.1, max_axis*1.1)

    surf = ax.plot_surface(x, y, z, linewidth=0, antialiased=False)
    temp = "" + str(v[0]) + "_" + str(v[1]) + "_" + str(v[2]) + "_" + str(v[3]) + "_" + str(v[4])
    plt.savefig("test" +temp+".pdf")




def print_all(v):
    print("Particle Values: a: {0:.4e} b: {1:.4e} c: {2:.4e} e1: {3:.4e} e2: {4:.4e}".format(v[0],v[1],v[2],v[3],v[4]))
    PMoI = moment_of_inertia(v)
    vol = volume(v)
    print("Volume: {0:.4e}  Prinical Momement of Inertia: [{1:.4e}, {2:.4e}, {3:.4e}]".format(vol, PMoI[0], PMoI[1], PMoI[2]))
    plot_and_save(v)

#v = [5.0e-3,5.0e-3,25.0e-6,0.2,0.2]



# v = [0.001,0.001,0.006,2.0/2.0,2.0/10.0]
# print_all(v)


#To double check  ellipsoid or spheres
# v[3] (e1) and v[4] (e2) should equal 1.0
# print(4.0/3.0 * math.pi * v[0]*v[1]*v[2])
# print((4.0 * math.pi / 15.0 * v[0] * v[1] * v[2] * (v[1]**2 + v[2]**2)))
# print((4.0 * math.pi / 15.0 * v[0] * v[1] * v[2] * (v[0]**2 + v[2]**2)))
# print((4.0 * math.pi / 15.0 * v[0] * v[1] * v[2] * (v[0]**2 + v[1]**2)))



#To double check cubes
# print(4.0/3.0 * math.pi * v[0]*v[1]*v[2])
# print((1.0 / 12.0 * v[0] * v[1] * v[2] * (v[1]**2 + v[2]**2)))
# print((1.0 / 12.0 * v[0] * v[1] * v[2] * (v[0]**2 + v[2]**2)))
# print((1.0 / 12.0 * v[0] * v[1] * v[2] * (v[0]**2 + v[1]**2)))


