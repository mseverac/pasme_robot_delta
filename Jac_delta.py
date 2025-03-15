import numpy as np
from MGD_delta import *
from MGI_delta import *

def Jac_delta(alpha, l1, l2, r,X=np.array([])):
    """calcule les jacobiennes tq J1*dX=J2*dQ"""
    phi1=2*np.pi/3
    phi2=4*np.pi/3
    c1=np.cos(phi1)
    s1=np.sin(phi1)
    c2=np.cos(phi2)
    s2=np.sin(phi2)
    R01 = np.matrix([[c1, -s1, 0], [s1, c1, 0], [0, 0, 1]])
    R02 = np.matrix([[c2, -s2, 0], [s2, c2, 0], [0, 0, 1]])
    R00 = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    R10 = np.matrix([[c1, s1, 0], [-s1, c1, 0], [0, 0, 1]])
    R20 = np.matrix([[c2, s2, 0], [-s2, c2, 0], [0, 0, 1]])
    R00 = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    pis2=np.pi/2
    c3=np.cos(pis2)
    s3=np.sin(pis2)
    E=np.matrix([[c3, 0, s3], [0, 1, 0], [-s3, 0, c3]])
    Rs=[R00,R01,R02]
    Rsi=[R00,R10,R20]
    Cs,us=delta_pts_C_u(alpha,l1,r)
    # print(f"alpha = {alpha}")
    if len(X)==0 : X=MGD_delta(alpha, l1, l2, r)
    # print(f"X = {X}")
    Ds=delta_pts_D(X,r)
    vs=[]

    Eqs=[]
    Equs=[]
    Equvl=[]
    for i in range(3):
        # print(f"Cs[{i}] = {Cs[i]}")
        # print(f"Ds[{i}] = {Ds[i]}")
        v=Ds[i]-Cs[i]
        v=(v/np.linalg.norm(v))[0]
        vs.append(v)

        Eq=np.array(Rs[i]@E)
        Eqs.append(Eq)

        Equs.append((Eq@Rsi[i]@us[i])[0])
        Equvl.append((np.dot(Equs[i],v))[0,0])


    J1=np.matrix([vs[0],vs[1],vs[2]])

    J2=np.matrix([[Equvl[0],0,0], [0, Equvl[1], 0], [0, 0, Equvl[2]]])

    return J1,J2



l1 = 0.260
l2 = 0.520
r = 0.634 * l1
r2 = 0.033
r1 = r + r2
