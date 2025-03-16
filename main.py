import numpy as np
from scipy.optimize import minimize
from MGD_delta import *
from MGI_delta import *
from Jac_delta import *
from condition import *
import matplotlib.pyplot as plt


Kp=1   #gain pour optimisation du poid
Kt=0.1 #gain pour optimisation de la taille


# Fonction objectif
def objective(x):
    l1,l2,r1,r2,a,b=x

    A=a**2-b**2
    V=A*(l1+l2)
    T=l1+l2
    return Kp*V+Kt*T
  # Correct the typo from 'somme' to 'sum'

# Bornes des variables
bounds = [(0.1, 10), (0.1, 10), (0.1, 10), (0.033, 10),(0.01,0.3),(0.0,0.1)]

constraints = [{'type': 'ineq', 'fun': lambda x, i=i: -condition(x)[i]} for i in range(len(condition(np.zeros(6))))]  # Adjust the range to the size of the condition array

# Options pour l'optimisation
options = {
    'maxiter': 5,       # Nombre maximal d'itérations (par défaut : 100)
    'iprint': 2,          # Niveau de détail des affichages (0 = silence, 1 = affichage de base, 2+ = détaillé)
    'disp': True,         # Afficher les messages de convergence (True/False)
}
# Point initial
x0 = np.array([0.6, 1.2, 0.634 * 0.6, 0.033,0.1,0.05])

# Optimisation
result = minimize(objective, x0, method='SLSQP', bounds=bounds, constraints=constraints,options=options)

# Résultat
print("Solution optimale :", result.x)
print("Valeur minimale de la fonction :", result.fun)
