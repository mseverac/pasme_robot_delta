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


"""print("b (radians):", b)
print("Cs:", Cs)
print("Ds:", Ds)
print("us:", us)
print("b (degres):", b_degre)
print("X recalcule:", X_calculated)"""


import matplotlib.pyplot as plt

# Génération des points sur l'axe z
z_values = np.linspace(-0.35, -0.75, 100)
condition_numbers = []

for z in z_values:
    X = np.array([0.0, 0.0, z])
    b = MGI_delta(X, l1, l2, r)
    J1,J2 = Jac_delta(b, l1, l2, r)
    J=np.linalg.inv(J1)@J2
    print(J)
    condition_number = np.linalg.cond(J)
    condition_numbers.append(condition_number)

# Création de la figure
plt.figure()
plt.plot(z_values, condition_numbers, label='Condition number of J')
plt.xlabel('z')
plt.ylabel('Condition number')
plt.title('Condition number of the Jacobian matrix J for points on the z-axis')
plt.legend()
plt.grid(True)
plt.show()

