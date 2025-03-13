import numpy as np
from MGD_delta import *
from MGI_delta import *
from Jac_delta import *

# Définition des constantes
l1 = 0.260
l2 = 0.520
r = 0.634 * l1
r2 = 0.033
r1 = r + r2

# Constantes supplémentaires
gammaMax = 40 * (np.pi / 180)
alphaBetaMin = 30 * (np.pi / 180)

D = 2 * l2 * np.sin(gammaMax)
H = l1 * np.sqrt((1 + l2 / l1) ** 2 - (r / l1 + 2 * ((l2 * np.sin(gammaMax)) / (l1 * np.sqrt(3))))) - l1 * np.sqrt(1 + (l2 / l1) ** 2 - (l2 / l1) * np.sqrt(3))
Zh = -l1 * np.sqrt(1 + (l2 / l1) ** 2 - 2 * (l2 / l1) * np.cos(alphaBetaMin))

# Test des fonctions
X = np.array([0.0, 0.0, -0.35])
b = MGI_delta(X, l1, l2, r)
b_degre = (180 / np.pi) * b
X_calculated = MGD_delta(b, l1, l2, r)
Cs,us=delta_pts_C_u(b,l1,r1)
Ds=delta_pts_D(X_calculated,r2)


print("b (radians):", b)
print("Cs:", Cs)
print("Ds:", Ds)
print("us:", us)
print("b (degres):", b_degre)
print("X recalcule:", X_calculated)


