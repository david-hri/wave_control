

from scipy import signal
from math import *
from data import *
import cmath
from scipy.integrate import quad
from matplotlib import pyplot as plt
from scipy.optimize import minimize_scalar,minimize
# on notera w pour omega
def g1(y,w,amplitude=1.0):  # gausienne pour modéliser un bruit localiser
    try:
        if -L <= y <= L:
            return amplitude * sin(w * y)
        else:
            raise ValueError("Error: y is outside the valid range [-L, L]")
    except ValueError as e:
        return str(e)


def g2(y,w, amplitude=1.0): # 1 si w in [20,20000], pour modéliser une salle de classe
    try:
        if -L<=y<=L:
            freq=w/(2*pi)
            if 20<=freq<=20000:
                return amplitude 
            else: return 0
        else: 
            raise ValueError("Error: y is outside the valid range [-L, L]")
    except ValueError as e:
        return str(e)

def g3(y,w, amplitude =1.0,w0=10000,sigma=1000,k_sirène=1): #sirène , sous la forme d'un signal carr;e
     #la sirène se trouve entre -k_sirène*L et k_sirène*L
    if -k_sirène*L<=y<=k_sirène*L:
        return amplitude/(sqrt(2*pi)*sigma) *exp(-(w-w0)**2/sigma**2)
    else:
        return 0
    




def compute_fourier_coefficient(g,k, w): #calcul le coefficient k de la transformée de fourrier de g 
    integrand = lambda y: g(y, w) * np.exp(-1j * k  * y)
    result, _ = quad(integrand, -L, L)
    ck = complex(result.real / (2 * L), result.imag / (2 * L))
    return ck





def calculate_e_k2(alpha, w, k,g,A=1,B=1): #calcul des ek avec la formule de l'énoncé
    lambda0=f_lambda0(k,w)
    lambda1=f_lambda1(k,w)
    gk=compute_fourier_coefficient(g,k,w)
    chi=f_chi(k, lambda0, lambda1, gk,alpha)
    gamma=f_gamma(k, alpha, lambda0,lambda1, gk)
    if k**2>=Ksi0/eta0*w**2:
        result=(A+B*k**2)*(1/(2*lambda0)*chi*chi.conjugate()*(1+np.exp(-2*lambda0*L))+gamma*gamma.conjugate()*(-1+np.exp(2*lambda0*L))+2*L*((chi*gamma.conjugate()).real))
        result+=B*lambda0/2*chi*chi.conjugate()*(1+np.exp(-2*lambda0*L))+gamma*gamma.conjugate()*(-1+np.exp(2*lambda0*L))-2*B*lambda0**2*L*((chi*gamma.conjugate()).real)
        return result
    else: 
        result=chi*gamma.conjugate()*(1+np.exp(-2*lambda0*L)).imag
        result/=lambda0
        result=(A+B*k**2)*complex(L*(chi*chi.conjugate()+gamma*gamma.conjugate()),result)
        result+=B*lambda0*lambda0.conjugate()*L*(chi*chi.conjugate()+gamma*gamma.conjugate())
        result+=complex(0,B*lambda0*(chi*gamma.conjugate()*(1+np.exp(-2*lambda0*L))).imag)
        r1=L*(chi*chi.conjugate()+gamma*gamma.conjugate())
        i1=chi*gamma.conjugate()*(1+np.exp(-2*lambda0*L))
        i1=i1.imag
        i3=B*lambda0*i1
        i1=1/lambda0*i1
        temp1=complex(r1,i1)
        temp1=(A+B*k**2)*temp1
        temp2=B*lambda0*lambda0.conjugate()*r1
        temp3=complex(0,i3)
        return temp1+temp2+temp3
        # return result


def e(alpha,w,N,g): #somme des ek
    result=0
    for i in range(-N,N+1):
        k=i*np.pi/L
        result+=calculate_e_k2(alpha,w,k,g)
    return result




def optimal_alpha(w,N,g): #utilisation du module scipy pour trouver le alpha optimal
    init=[w/10,w/10]
    def optimize_e(xy):
        x,y=xy
        return abs(e(x+1j*y,w,N,g))
    result = minimize(optimize_e,init,method='L-BFGS-B')
    return result






def minimise(g,W): #fonction principale qui donne alpha optimal pour les différents w
    initial_alpha=[1,1.5]
    alpha_real=[]
    alpha_img=[]
    for w in W:
        print(w)
        def objective_function(xy):
            x,y=xy
            res=e(x+1j*y,w,20,g)
            return np.abs(res)
        result=minimize(objective_function,initial_alpha,method='L-BFGS-B')
        alpha_real.append(result.x[0])
        alpha_img.append(result.x[1])
    

    return alpha_real,alpha_img







