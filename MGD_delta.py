import numpy as np

def MGD_delta(alpha, l1, l2, r):
    # Définition des constantes
    phi = np.array([0, (2 * np.pi) / 3, (4 * np.pi) / 3])
    
    # Calcul de constantes utiles
    D = -l2**2 + l1**2 + r**2 + 2 * r * l1 * np.cos(alpha)
    D1, D2, D3 = D
    
    E = 2 * (r + l1 * np.cos(alpha)) * np.cos(phi)
    E1, E2, E3 = E
    
    F = 2 * (r + l1 * np.cos(alpha)) * np.sin(phi)
    F1, F2, F3 = F
    
    G = -2 * l1 * np.sin(alpha)
    G1, G2, G3 = G
    
    H1 = -(E3 - E1) * (G2 - G1) + (G3 - G1) * (E2 - E1)
    H2 = -(F3 - F1) * (E2 - E1) + (E3 - E1) * (F2 - F1)
    H3 = -(D3 - D1) * (E2 - E1) + (E3 - E1) * (D2 - D1)
    H4 = (D3 - D1) * (F2 - F1) - (F3 - F1) * (D1 - D2)
    H5 = -(G3 - G1) * (F2 - F1) - (F3 - F1) * (G1 - G2)
    
    L = 1 + (H5**2 + H1**2) / (H2**2)
    M = -(2 * (H5 * H4 + H1 * H3)) / (H2**2) + (E1 * H5 + F1 * H1) / H2 + G1
    N = D1 + (H4**2 + H3**2) / (H2**2) - (E1 * H4 + F1 * H3) / H2
    
    # Calcul des positions x, y, z
    z = (M - np.sqrt(M**2 - 4 * L * N)) / (2 * L)
    x = (H5 / H2) * z + H4 / H2
    y = (H1 / H2) * z + H3 / H2
    
    return np.array([x, y, z])



def delta_pts_C_u(alphas,l1,r1):
    phi1=2*np.pi/3
    phi2=4*np.pi/3
    c1=np.cos(phi1)
    s1=np.sin(phi1)
    c2=np.cos(phi2)
    s2=np.sin(phi2)
    R01 = np.matrix([[c1, -s1, 0], [s1, c1, 0], [0, 0, 1]])
    R02 = np.matrix([[c2, -s2, 0], [s2, c2, 0], [0, 0, 1]])
    R00 = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    Rs=[R00,R01,R02]
    Cs=[]
    us=[]
    B=np.array([r1,0,0])
    for i,alpha in enumerate(alphas) :
        u=np.array([l1*np.cos(alpha),0,-l1*np.sin(alpha)])
        C=B+u
        Cs.append(np.array((Rs[i]@C)[0]))
        u=np.array((Rs[i]@u))[0]
        u=(u/np.linalg.norm(u))
        us.append(u)
    return Cs,us

def delta_pts_D(MGD,r2):
    phi1=2*np.pi/3
    phi2=4*np.pi/3
    c1=np.cos(phi1)
    s1=np.sin(phi1)
    c2=np.cos(phi2)
    s2=np.sin(phi2)
    R01 = np.matrix([[c1, -s1, 0], [s1, c1, 0], [0, 0, 1]])
    R02 = np.matrix([[c2, -s2, 0], [s2, c2, 0], [0, 0, 1]])
    R00 = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    Rs=[R00,R01,R02]
    Ds=[]
    for R in Rs :
        Ds.append(np.array((R@(MGD+np.array([r2,0,0])))[0]))

    return Ds