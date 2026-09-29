import numpy as np
kB=1.380649e-23
T=273
J=0
H=1
mu=1
def Random_Grid_Generator_TwoD(N):
    return np.random.choice([-1,1],size=(N,N))

def Metropolis(dE):
    if dE >= 0:
        return np.exp(-dE/(kB*T))
    else:
        return 1
def Delta_Energy(grid,i,j):
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
    return 2*grid[i,j]*(J*neigbourSum+H*mu)

def Monte_Carlo_sim(N,nSim):
    lattice = Random_Grid_Generator_TwoD(N)
    Nx,Ny =lattice.shape
    randomI= np.random.randint(Nx,size=nSim)
    randomJ= np.random.randint(Ny,size=nSim)
    x = np.random.rand(nSim)
    for i in range(nSim):
        dE = Delta_Energy(lattice,randomI[i],randomJ[i])
        omega=Metropolis(dE)
        lattice[randomI[i],randomJ[i]]*= -1 if x[i]< omega else 1
    return lattice



