import numpy as np

def MGI_delta(X, l1, l2, r):
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
    t1 = (-A + np.sqrt(A**2 + B1**2 - C1**2)) / (B1 + C1)
    t2 = (-A + np.sqrt(A**2 + B2**2 - C2**2)) / (B2 + C2)
    t3 = (-A + np.sqrt(A**2 + B3**2 - C3**2)) / (B3 + C3)
    t = np.array([t1, t2, t3])
    
    alpha = 2 * np.arctan(t)
    
    return alpha

 