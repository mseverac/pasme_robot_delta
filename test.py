import numpy as np
from MGD_delta import *
from MGI_delta import *
from Jac_delta import *
from condition import *
import matplotlib.pyplot as plt


"il faut l1>0.37 si l2=2*l1 pour que le ws soit assez grand"

"il faut l1>0.6 si l2=2*l1 pour le conditionnement de la jacobienne"


"sans raideur sol optimale : [0.34790233 0.69230797 0.10000007 0.033     ]"


# Définition des constantes
l1 = 0.6
l2 = 1.2
r = 0.634 * l1
r2 = 0.033
r1 = r + r2
# Dictionnaire contenant les constantes du robot
robot_constants = {
    'l1': l1,
    'l2': l2,
    'r1': r1,
    'r2': r2
}

# Constantes supplémentaires
gammaMax = 40 * (np.pi / 180)
alphaBetaMin = 30 * (np.pi / 180)

D = 2 * l2 * np.sin(gammaMax)
H = l1 * np.sqrt((1 + l2 / l1) ** 2 - (r / l1 + 2 * ((l2 * np.sin(gammaMax)) / (l1 * np.sqrt(3))))) - l1 * np.sqrt(1 + (l2 / l1) ** 2 - (l2 / l1) * np.sqrt(3))
Zh = -l1 * np.sqrt(1 + (l2 / l1) ** 2 - 2 * (l2 / l1) * np.cos(alphaBetaMin))









































"""

C=condition([0.35704558, 0.66003772 ,0.1    ,    0.033     ],gammaMax,alphaBetaMin)
print(C)
print(max(C))
print(len(C))
print("condition faite")"""



























"""print("b (radians):", b)
print("Cs:", Cs)
print("Ds:", Ds)
print("us:", us)
print("b (degres):", b_degre)
print("X recalcule:", X_calculated)"""



"""X = np.array([0,0,-0.5])

# Définir la plage de valeurs pour l1
l1_values = np.linspace(0.1, 2, 100)
l2_values = l1_values * 2

# Calculer les valeurs de C de condition ws pour chaque l1
C_values = []
for l1, l2 in zip(l1_values, l2_values):
    r = 0.634 * l1
    r2 = 0.033
    r1 = r + r2
    b = MGI_delta(X, l1, l2, r)
    C, H, D, Zh = condition_ws(l1, l2, r1, r2, gammaMax, alphaBetaMin)
    C_values.append(C)

# Tracer la figure
C1_values = [C[0] for C in C_values]
C2_values = [C[1] for C in C_values]

plt.figure()
plt.plot(l1_values, C1_values, label='C1 de condition ws')
plt.plot(l1_values, C2_values, label='C2 de condition ws')
plt.xlabel('l1')
plt.ylabel('C de condition ws')
plt.title('C1 et C2 de condition ws en fonction de l1')
plt.legend()
plt.grid(True)
plt.show()"""

