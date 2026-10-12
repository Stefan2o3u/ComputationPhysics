import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from numba import njit
kB=1.380649e-23

def Random_Grid_Generator_TwoD(N):
    return np.random.choice([-1,1],size=(N,N))
def Up_Grid_Generator_TwoD(N):
    return np.ones((N,N))

@njit(cache=True)
def Metropolis(dE):
    if dE >= 0:
        return np.exp(-dE)
    else:
        return 1.0
@njit(cache=True)
def Delta_Energy(grid,i,j,Jdimless,Hdimless):
    """
    Algorithm for finding the change in energy of a single spin. for the neigbours I implemented
    periodic BC
    :param grid:
    :param index:
    :return:
    """
    Nx,Ny=grid.shape
    neigbourSum= (grid[(i+1) % Nx,j]+grid[(i-1) % Nx,j]
                  +grid[i,(j+1) % Ny]+grid[i,(j-1) % Ny])
    return 2*grid[i,j]*(Jdimless*neigbourSum+Hdimless)
@njit(cache=True)
def Sweep(lattice,J,H,nSweeps):
    N = lattice.shape[0]

    Mag = np.zeros(nSweeps)
    for i in range(nSweeps):
        I = np.random.randint(0,N, N * N)
        Jj = np.random.randint(0,N, N * N)
        x = np.random.rand(N * N)
        for k in range(N * N):
            dE = Delta_Energy(lattice, I[k], Jj[k], J, H)
            if x[k] < Metropolis(dE):
                lattice[I[k], Jj[k]] *= -1
        Mag[i] = lattice.sum()
    return Mag

def Monte_Carlo_sim(N,nSweep,J,H,Grid_Generator):
    lattice = Grid_Generator(N)
    Mag =Sweep(lattice,J,H,nSweep)

    im =np.full((N, N), False)
    im[lattice==-1]=True

    #plt.figure(figsize=(16,16))
    #plt.imshow(im)
    #plt.show()

    return lattice, Mag
#Monte_Carlo_sim(250,10000000,1,0)
def Average_Magnetisation(N,Hdimless):
    return np.tanh(Hdimless)
def check_42(N,H):
    lattice, Mag = Monte_Carlo_sim(N,500000,0,H,Random_Grid_Generator_TwoD)
    analytical = Average_Magnetisation(N,H)
    return Mag[-1]/(N**2), analytical
def plot_different_H(H):
    """" Asignment C"""
    analytical=np.zeros(len(H))
    avg = np.zeros(len(H))
    for i in range(len(H)):
        avg[i],analytical[i]=check_42(200,H[i])
    fig, ax = plt.subplots(figsize=(16,8))
    ax.plot(H,avg,label="Monte_Carlo")
    ax.plot(H,analytical,label="Analytical")
    ax.set(xlabel="muH/kBT", ylabel="<M>")
    # title=f"Convergence of {namep}-MovingMesh")
    ax.tick_params(axis="both", which="both", direction="in")
    ax.legend(loc="upper right")
    ax.grid()
    plt.tight_layout()
    plt.show()
#plot_different_H(np.linspace(0,8,50))


def Equibritaion_time(J,N,H,Generator,nSim):
    lattice,Mag=Monte_Carlo_sim(N,nSim,J,H,Generator)
    fig, ax = plt.subplots(figsize=(16, 8))
    ax.semilogx(np.linspace(1,nSim,nSim), Mag, label="Monte_Carlo")

    ax.set(xlabel="MC Sweeps", ylabel="<M>")
    # title=f"Convergence of {namep}-MovingMesh")
    ax.tick_params(axis="both", which="both", direction="in")
    ax.legend(loc="upper right")
    ax.grid()
    plt.tight_layout()
    plt.show()
    return None
#Equibritaion_time(0.1,32,0,Random_Grid_Generator_TwoD,10000)

def Onsager(J,N):
    M= np.zeros(len(J))
    M[J>1/2.269]=(1-np.sinh(2*J[J>1/2.269])**(-4))**(1/8)
    return M
def check_43(N,J):
    lattice, Mag = Monte_Carlo_sim(N,22000,J,0,Up_Grid_Generator_TwoD)
    return np.mean(Mag[2000:])/(N**2)
def Part_E():
    N=32
    J= np.linspace(0.1,2,1000)
    analytical = Onsager(J,N)
    avg = np.zeros(len(J))
    for i in range(len(J)):
        avg[i] = check_43(N, J[i])
    fig, ax = plt.subplots(figsize=(16, 8))
    ax.plot(J, avg, label="Monte_Carlo")
    ax.plot(J, analytical, label="Analytical")
    ax.set(xlabel="muH/kBT", ylabel="<M>")
    # title=f"Convergence of {namep}-MovingMesh")
    ax.tick_params(axis="both", which="both", direction="in")
    ax.legend(loc="upper right")
    ax.grid()
    plt.tight_layout()
    plt.show()
    return None
Part_E()
def Magnetic_Susceptibility(N,M):
    print(M)
    return 1/N*np.mean((M-np.mean(M))**2)
def Part_F():
    N= np.array([4, 6, 8, 10, 14, 20, 28, 40])
    J= 0.1
    J2=0.43
    xi = np.zeros((len(N),2))
    for i in range(len(N)):
        lat,M = Monte_Carlo_sim(N[i],21000,J,0,Up_Grid_Generator_TwoD)
        xi[i,0] = Magnetic_Susceptibility(N[i]**2,M[1000:])
        lat, M2 = Monte_Carlo_sim(N[i], 21000, J2, 0, Up_Grid_Generator_TwoD)
        xi[i, 1] = Magnetic_Susceptibility(N[i] ** 2, M2[1000:])
    fig, ax = plt.subplots(figsize=(16, 8))
    ax.loglog(N, xi[:,0], label="Monte_Carlo")
    ax.loglog(N, xi[:,1], label="Monte_Carlo")

    ax.set(xlabel="Size", ylabel="xi")
    # title=f"Convergence of {namep}-MovingMesh")
    ax.tick_params(axis="both", which="both", direction="in")
    ax.legend(loc="upper right")
    ax.grid()
    plt.tight_layout()
    plt.show()


    return None
#Part_F()