import numpy as np
import matplotlib.pyplot as plt
def P_function(A,phi):
    P = 1/(2*np.pi)*(1+A*np.cos(phi))
    return P

def rejection(A):
    y=np.random.uniform(0,1/np.pi,10000)
    phi = np.random.uniform(0,2*np.pi, 10000)
    Pvaues = P_function(A,phi)
    accepted_phi= phi[y <=Pvaues]

    return accepted_phi
#rejection(1)
def potter(A):
    phi = np.linspace(0,2*np.pi,1000)
    P = P_function(0.2,phi)
    plt.plot(P,phi)
    plt.show()