from functions import *
import matplotlib.pyplot as plt
X=np.linspace(-L,L,1000)
Y=[]
for i in range(len(X)):
    print(i)
    Y.append(abs(e(X[i],1000,20,g3)))
plt.plot(X,Y)
plt.show()




import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_3d_graph(X, Y, Z, title='Graphe 3D', xlabel='X', ylabel='Y', zlabel='Z'):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)

    ax.plot_surface(X, Y, Z, cmap='viridis')

    plt.show()

# Exemple d'utilisation :
import numpy as np

# Générer des données pour X, Y et Z (remplacez cela par vos données)
X, Y = np.meshgrid(np.linspace(-L, L, 100), np.linspace(-L, L, 100))
Z = np.zeros(X.shape)
for i in range(X.shape[0]):
    for k in range(X.shape[1]):
        # Calculer Z = e_alpha(X + iY) (dans cet exemple, Z = sin(sqrt(X^2 + Y^2))
        Z[i, k] = abs(e(complex(X[i,k],Y[i,k]),10000,10,g3))


#plot_3d_graph(X, Y, Z, 'Graphe 3D de sin(sqrt(X^2 + Y^2))', 'X', 'Y', 'Z')