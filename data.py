import numpy as np
from math import *
L=10

#dictionnaire avec les valeurs de phi sigma et alpha_h pour les différents matériau
dictionnaire = {
    "foam": [0.99, 14000, 1.02],
    "pine": [0.54, 35540, 1.22],
    "oak": [0.67, 197525, 1.36]
}


#toutes ces constantes sont données dans l'énoncé
phi,sigma,alpha_h=dictionnaire["oak"]
rho0=1.29

gamma_p=7/5
c0=340
Ksi0=1/c0/c0
Ksi1=Ksi0*gamma_p*phi
a0=0
a1=sigma*phi**2*gamma_p/(c0**2*rho0*alpha_h)
eta0=1
eta1=phi/alpha_h

def f_lambda0(k, omega): #calcul de lambda 0 pour k et w donné
    try:
        if k**2 >= Ksi0 * eta0 * omega**2:
            return np.sqrt(k**2 - Ksi0 / eta0 * omega**2)
        else:
            return 1j * np.sqrt(Ksi0 /eta0 * omega**2 - k**2)
    except:
        print(k,omega,"erreur")
    

def f_lambda1(k,w): #calcul de lambda 0 pour k et w donné
    k2=k**2 
    w2=Ksi1/eta1*w**2
    diff=k2-w2
    var1=(a1/eta1*w)**2 
    var2=sqrt(diff**2+var1)
    real=sqrt(diff+var2)
    real=1/sqrt(2)*real
    im=sqrt(-diff+var2)
    im=-1/sqrt(2)*im
    return complex(real,im)


# Fonction f(x)
def f( lambda0,x):
    return (lambda0 * eta0 - x) * np.exp(-lambda0 * L) + (lambda0 * eta0 + x) * np.exp(lambda0 * L)

# Fonction χ(k, α)
def f_chi(k, lambda0, lambda1, gk,alpha):
    result = gk * ((lambda0 * eta0 - lambda1 * eta1)/f(lambda0,lambda1*eta1)- (lambda0 * eta0 - alpha)/f(lambda0,alpha))
    return result

# Fonction γ(k, α)
def f_gamma(k, alpha, lambda0,lambda1, gk):
    numerator = gk * ((lambda0 * eta0 + lambda1 * eta1)/f(lambda0,lambda1*eta1)- (lambda0 * eta0 + alpha)/f(lambda0,alpha))
    return numerator


