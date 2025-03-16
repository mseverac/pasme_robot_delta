import numpy as np
from MGD_delta import *
from MGI_delta import *
from Jac_delta import *

gammaMax = 40 * (np.pi / 180)
alphaBetaMin = 30 * (np.pi / 180)

def condition_ws(l1,l2,r1,r2,gammaMax = 40 * (np.pi / 180),alphaBetaMin = 30 * (np.pi / 180)):
    """condition que le ws soit assez grand, C,H,D,Zh"""

    r=r1-r2
    D = 2 * l2 * np.sin(gammaMax)
    H = l1 * np.sqrt((1 + l2 / l1) ** 2 - (r / l1 + 2 * ((l2 * np.sin(gammaMax)) / (l1 * np.sqrt(3))))) - l1 * np.sqrt(1 + (l2 / l1) ** 2 - (l2 / l1) * np.sqrt(3))
    Zh = -l1 * np.sqrt(1 + (l2 / l1) ** 2 - 2 * (l2 / l1) * np.cos(alphaBetaMin))

    C=np.array([0.51-H,0.6-D/np.sqrt(2)])

    return C,H,D,Zh


def discrétisation(centre,pas):
    """renvoie une list de points 3D pour le workspace de coté 500mm autour du centre espacés de pas m"""
    n=int(0.5/pas)
    step=0.5/n
    points=[]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                points.append(np.array([centre[0]+i*step-0.25,centre[1]-0.25+j*step,centre[2]-0.25+k*step]))
    return points


def condition_raideur(Ktet,J):
    Cx=J@np.linalg.inv(Ktet)@J.T
    normCx=np.linalg.norm(Cx,2)
    dpmax=normCx*10
    print(f"Ktet = {Ktet}")
    print(f"normCx = {normCx}")
    print(f"dpmax = {dpmax}")
    print(f"Cx = {Cx}")
    return (0.0001-dpmax)*1000




def condition(x,gammaMax = 40 * (np.pi / 180),alphaBetaMin = 30 * (np.pi / 180)):
    """renvoie True si les conditions sont vérifiées, False sinon"""
    l1,l2,r1,r2,a,b=x
    Cgeo=np.array([b-a+0.01])*1000

    I=(pow(a,4)-pow(b,4))/12
    E=70e9
    k=3*E*I/(l1**3)


    Ktet=np.diag([k,k,k])

    Cws,H,D,Zh=condition_ws(l1,l2,r1,r2)
    centre = np.array([0,0,Zh-0.25])
    pts=discrétisation(centre,0.15)
    C_cond=[]
    C_raideur=[]
    if Cws[0]>0 or Cws[1]>0:
        print("ws pas assez grand")
        print(f"l1 = {l1}")
        print(f"l2 = {l2}")
        print(f"sin gammaMax = {np.sin(gammaMax)}")
        print(f"Zh = {Zh}")
        print(f"H = {H}")
        print(f"D = {D}")
        print(len(np.concatenate((Cws, np.ones(2*len(pts))))))
        return np.concatenate((Cws, np.ones(len(pts))))
    for pt in pts:
        
        alpha=MGI_delta(pt,l1,l2,r1)
        J1,J2=Jac_delta(alpha,l1,l2,r1,X=pt)
        if np.linalg.det(J1)==0:   
            print("erreur de singularité sur J1")
            print(f"pt = {pt}")
            print(f"l1 = {l1}")
            print(f"l2 = {l2}")
            print(f"r1 = {r1}")
            print(f"r2 = {r2}")
            print(f"gammaMax = {gammaMax}")
            print(f"alphaBetaMin = {alphaBetaMin}")
            print(f"J1 = {J1}")
            print(f"J2 = {J2}")
            print(len(Cws))
            return Cws
        if np.linalg.det(J2)==0:
            print("erreur de singularité sur J2")
            print(f"pt = {pt}")
            print(f"l1 = {l1}")
            print(f"l2 = {l2}")
            print(f"r1 = {r1}")
            print(f"r2 = {r2}")
            print(f"gammaMax = {gammaMax}")
            print(f"alphaBetaMin = {alphaBetaMin}")
            print(f"J1 = {J1}")
            print(f"J2 = {J2}")
            print(len(Cws))
            return Cws
        J=np.linalg.inv(J1)@J2
        C_raideur.append(condition_raideur(Ktet,J))
       
        cond=np.linalg.cond(J)
        if cond==np.inf:
            C_cond.append(0.1)
        else :
            C_cond.append(0.1-1/cond)

    C_cond=np.array(C_cond)
    C_raideur=np.array(C_raideur)
    C = np.concatenate((Cgeo,Cws, C_cond,C_raideur))
    return C


    