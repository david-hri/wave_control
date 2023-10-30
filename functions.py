

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
    


def fourrier_k(g,w,k): #calcul le coefficient k de la transformée de fourrier de g 
    def gprime(y):
        return g(y,w) * np.exp(-complex(0,1) * np.pi * k * y /L) #le calcul du coefficient de fourrier nous fait calculer l'intégrale de cette fonction entre -L et L
    result, error = quad(gprime, -L, L)
    return result

def compute_fourier_coefficient(g,k, w):
    # Define the integrand for the Fourier coefficient
    integrand = lambda y: g(y, w) * np.exp(-1j * k  * y)

    # Perform the integration over the range [-L, L]
    result, _ = quad(integrand, -L, L)

    # Divide by 2*L to compute the coefficient
    ck = complex(result.real / (2 * L), result.imag / (2 * L))
    return ck


def calculate_e_k(alpha, w, k,g,A=1,B=1):
    lambda0=f_lambda0(k,w)
    lambda1=f_lambda1(k,w)
    gk=compute_fourier_coefficient(g,k,w)
    chi=f_chi(k, lambda0, lambda1, gk,alpha)
    gamma=f_gamma(k, alpha, lambda0,lambda1, gk)
    if k**2>=Ksi0/eta0*w**2:
        temp1=chi*chi.conjugate()*(1+np.exp(-2*lambda0*L))+gamma*gamma.conjugate()*(-1+np.exp(2*lambda0*L))
        temp2=temp1
        prod=chi*gamma.conjugate()
        temp1=1/(2*lambda0)*temp1+2*L*prod.real
        temp1=(A+B*k**2)*temp1
        temp2=B*lambda0/2*temp2-2*B*lambda0**2*L*prod.real
        print("ici")
        
    if k**2 >= Ksi0 / eta0 * w**2:
        result = (A + B * k**2) * (L * (abs(chi)**2 + abs(gamma)**2) +
                 (1j / lambda0) * (chi.conjugate() * gamma * (1 - cmath.exp(-2 * lambda0 * L))))
        result += B * L * abs(lambda0)**2 * (abs(chi)**2 + abs(gamma)**2) + 1j * B * lambda0 * (chi.conjugate() * gamma * (1 - cmath.exp(-2 * lambda0 * L)))
    
    
    else:
        result = (A + B * k**2) * (1 / (2 * lambda0) * ((abs(chi)**2 * (1 - cmath.exp(-2 * lambda0 * L))) + (abs(gamma)**2 * (cmath.exp(2 * lambda0 * L) - 1))) + 2 * L * (chi.real * gamma.real))
        result += B * (lambda0 / 2) * ((abs(chi)**2 * (1 - cmath.exp(-2 * lambda0 * L)) + abs(gamma)**2 * (cmath.exp(2 * lambda0 * L) - 1)) - 2 * B * lambda0**2 * L * (chi.real * gamma.real))
    return result



def calculate_e_k2(alpha, w, k,g,A=1,B=1):
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


def e(alpha,w,N,g):
    result=0
    for i in range(-N,N+1):
        k=i*np.pi/L
        result+=calculate_e_k2(alpha,w,k,g)
    return result




def optimal_alpha(w,N,g):
    init=[w/10,w/10]
    def optimize_e(xy):
        x,y=xy
        return abs(e(x+1j*y,w,N,g))
    result = minimize(optimize_e,init,method='L-BFGS-B')
    return result






def minimise(g,W): #modify the function with the correct source, do not forget to modify the labels as well
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
    
    #plt.plot(w_range, alpha_real, color='red', label=f'partie Reelle de alpha pour {g.__name__} ')  
    return alpha_real,alpha_img








