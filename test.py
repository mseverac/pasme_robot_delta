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










































x=[4.41053309e-01 ,1.02938724e+00, 3.74043341e-01, 3.30000000e-02,1.00000000e-02 ,3.49751378e-17]
C=condition(x)
print(C)
print(C[0])
print(C[1])
print(C[2])
print(C[3])
print(C[4])
print(C[5])
print(C[6])
print(C[7])
print(C[8])
print(C[9])
print(C[10])
print(C[11])
print(C[12])
print(C[13])
print(C[14])
print(C[15])
print(C[16])
print(C[17])
print(C[18])
print(C[19])
print(C[20])
print(C[21])
print(max(C))
print("Index of max value in C:", np.argmax(C))
print(C[np.argmax(C)])
print(C[np.argmax(C)+1])
print(C[np.argmax(C)+2])

print("condition faite")



























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

