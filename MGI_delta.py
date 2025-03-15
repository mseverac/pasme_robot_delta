import numpy as np


def condition_ws(l1,l2,r,gammaMax = 40 * (np.pi / 180),alphaBetaMin = 30 * (np.pi / 180)):
    """condition que le ws soit assez grand, C,H,D,Zh"""

    
    D = 2 * l2 * np.sin(gammaMax)
    H = l1 * np.sqrt((1 + l2 / l1) ** 2 - (r / l1 + 2 * ((l2 * np.sin(gammaMax)) / (l1 * np.sqrt(3))))) - l1 * np.sqrt(1 + (l2 / l1) ** 2 - (l2 / l1) * np.sqrt(3))
    Zh = -l1 * np.sqrt(1 + (l2 / l1) ** 2 - 2 * (l2 / l1) * np.cos(alphaBetaMin))

    C=np.array([0.5-H,0.5-D/np.sqrt(2)])

    return C,H,D,Zh


def MGI_delta(X, l1, l2, r):
    """Calcule les angles alpha pour un point donné"""
    cws,H,D,Zh = condition_ws(l1, l2, r, )
    R=np.sqrt(X[0]**2+X[1]**2)
    if X[2]>Zh or X[2]<Zh-H or R>D/2:
        print("erreur de ws")
        print(f"Zh = {Zh}")
        print(f"H = {H}")
        print(f"D = {D}")
        
        raise ValueError(f"\033[91mPoint hors workspace {X}\033[0m")
        return [None, None, None]
    # Définition des constantes
    phi = np.array([0, (2 * np.pi) / 3, (4 * np.pi) / 3])
    
    # On récupère les coordonnées
    x, y, z = X
    
    # Calcul de constantes utiles
    A = -2 * l1 * z
    B = 2 * l1 * (r - x * np.cos(phi) - y * np.sin(phi))
    B1, B2, B3 = B
    
    C = l2**2 - l1**2 - r**2 - x**2 - y**2 - z**2 + 2 * r * (x * np.cos(phi) + y * np.sin(phi))
    C1, C2, C3 = C
    
    # Calcul des alphas
    if A**2 + B1**2 - C1**2 < 0 or A**2 + B2**2 - C2**2 < 0 or A**2 + B3**2 - C3**2 < 0:
        raise ValueError(f"\033[91mPoint hors workspace {X}\033[0m")
        return [None, None, None]
    t1 = (-A + np.sqrt(A**2 + B1**2 - C1**2)) / (B1 + C1)
    t2 = (-A + np.sqrt(A**2 + B2**2 - C2**2)) / (B2 + C2)
    t3 = (-A + np.sqrt(A**2 + B3**2 - C3**2)) / (B3 + C3)
    t = np.array([t1, t2, t3])
    
    alpha = 2 * np.arctan(t)
    
    return alpha

 