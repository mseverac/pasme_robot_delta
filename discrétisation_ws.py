

def discrétisation(centre,pas):
    """renvoie une list de points 3D pour le workspace de coté 500mm autour du centre espacés de pas mm"""
    n=int(500/pas)
    step=500/n
    points=[]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                points.append([centre[0]+i*step,centre[1]+j*step,centre[2]+k*step])
    return points

