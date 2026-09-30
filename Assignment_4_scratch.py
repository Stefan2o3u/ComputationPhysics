import numpy as np
import matplotlib.pyplot as plt
kB=1.380649e-23

def Random_Grid_Generator_TwoD(N):
    return np.random.choice([-1,1],size=(N,N))

def Metropolis(dE):
    if dE >= 0:
        return np.exp(-dE)
    else:
        return 1
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

def Monte_Carlo_sim(N,nSim,J,H):
    lattice = Random_Grid_Generator_TwoD(N)
    Nx,Ny =lattice.shape
    randomI= np.random.randint(Nx,size=nSim)
    randomJ= np.random.randint(Ny,size=nSim)
    x = np.random.rand(nSim)
    print(lattice)
    for i in range(nSim):
        dE = Delta_Energy(lattice,randomI[i],randomJ[i],J,H)
        omega=Metropolis(dE)

        lattice[randomI[i],randomJ[i]]*= -1 if x[i]< omega else 1
    im =np.full((N, N), False)
    im[lattice==-1]=True
    print(im)
    plt.figure(figsize=(16,16))
    plt.imshow(im)
    plt.show()
    return lattice
Monte_Carlo_sim(250,10000000,1,0)
def Average_Magnetisation(N):
    return mu*N*np.tanh(mu*H/(kB*T))
def check_42(N):
    lattice = Monte_Carlo_sim(N,100000)
    print(mu*np.average(lattice))
    analytical = Average_Magnetisation(N)
    print(analytical)
#check_42(200)
