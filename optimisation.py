import numpy as np
from scipy.optimize import minimize
from main import *


# Fonction objectif
def objective(x):
    return np.sum(x)  # Correct the typo from 'somme' to 'sum'

# Bornes des variables
bounds = [(0.1, 10), (0.1, 10), (0.1, 10), (0.033, 10)]

constraints = [{'type': 'ineq', 'fun': lambda x, i=i: -condition(x)[i]} for i in range(len(condition(np.zeros(4))))]  # Adjust the range to the size of the condition array

# Options pour l'optimisation
options = {
    'maxiter': 10,       # Nombre maximal d'itérations (par défaut : 100)
    'iprint': 2,          # Niveau de détail des affichages (0 = silence, 1 = affichage de base, 2+ = détaillé)
    'disp': True,         # Afficher les messages de convergence (True/False)
}
# Point initial
x0 = np.array([0.6, 1.2, 0.634 * 0.6, 0.033])

# Optimisation
result = minimize(objective, x0, method='SLSQP', bounds=bounds, constraints=constraints,options=options)

# Résultat
print("Solution optimale :", result.x)
print("Valeur minimale de la fonction :", result.fun)
