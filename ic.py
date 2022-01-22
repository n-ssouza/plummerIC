import os
import numpy as np
from configparser import ConfigParser


# Unit System: 
# [L] = kpc 
# [V] = km/s
# [M] = 1e10 Msun
# [T] = 0.98 Gyr
# [G] = kpc (km/s)^2 (1e10 Msun)^{-1}


def Rho(M, a, r):
    return 0.75*(M/(np.pi*(np.power(a, 3)))) * np.power((1 + (np.square(r) / np.square(a))), -2.5)

def Phi(M, a, r):
    G = 43007.1    # [L] = kpc ; [M] = 1e10 Msun ; [T] approx 0.98 Gyr ; [V] = km/s
    return -(G*M / np.sqrt(np.square(r)+np.square(a)))

def CumulativeMass(M, a, r):
    return M * np.power(1 + (np.square(a/r)), -1.5)

def InverseCumulativeMass(M, a):
    return a*(np.power((np.power(M, -2/3) - 1), -0.5))

def UniSphDistr(R, N):
    #Spherical coordinates: cos(phi) is randomly sampled in order for the sphere to be isotropic
    theta = 2*np.pi*np.random.rand(N)    #  0 < theta < 2pi
    cosphi = 2*np.random.rand(N) - 1     # -1 < cos(phi) < 1
    phi = np.arccos(cosphi)

    x = R * np.cos(theta) * np.sin(phi)
    y = R * np.sin(theta) * np.sin(phi)
    z = R * np.cos(phi)
    
    return x,y,z

def EscapeVel(M, a, R):
    return np.sqrt(2 * (-Phi(M, a, R)))

def VelDistribution(q):
    return np.square(q)*np.power(1 - np.square(q), 3.5)


def main():
    init()

    Positions, Rmags = set_positions()
    Velocities = set_velocities(Rmags)
    Masses = set_masses()

    write_file(Positions, Velocities, Masses)


def init(): 
    global Mt, a, N, Rtrunc, filename, filetype

    config = ConfigParser()
    config.read('plummer_params.ini')

    Mt = config.getfloat('plummer', 'Mt')          # 1e10 Msun
    a = config.getfloat('plummer', 'a')            # kpc
    Rtrunc = config.getfloat('plummer', 'Rtrunc')  # kpc
    N = config.getint('plummer', 'N')              # Number of particles
    
    filetype = config.get('global', 'filetype')
    filename = config.get('global', 'filename')


#---------------Positions-----------------#

#Cummulative mass / Total mass  -->  numbers between 0 and 1 

def set_positions():
    R = np.zeros(N)
    i = 0
    while (i < N):
        Mcumm = np.random.rand()
        radius = InverseCumulativeMass(Mcumm, a)
        if (radius < Rtrunc):
            R[i] = radius
            i = i + 1 

    return UniSphDistr(R, N), R 


#---------------Velocities----------------#

#von Neumann Rejection: randomly choose valid velocity magnitudes for each particle (given their radii) 

def set_velocities(R):
    validVelMag = np.zeros(N)
    i = 0
    while (i < N):
        #For each particle: define Q = V / Vescape, this way 0 < Q < 1
        Q = np.random.rand()
        DistrQ = 0.1 * np.random.rand()        # 0.1 is the maximum of VelDistributio function

        if (DistrQ < VelDistribution(Q)):      #A valid (Q, DistrQ) pair is found:
            validVelMag[i] = Q * EscapeVel(Mt, a, R[i])
            i = i + 1

    return UniSphDistr(validVelMag, N)


#-----------------Masses------------------#

def set_masses():
    return np.ones(N) * (CumulativeMass(Mt, a, Rtrunc) / N)


#------------Saving-Text-File-------------#

def write_file(Positions, Velocities, Masses):
    
    savePath = './InitialConditions/'
    
    if not (os.path.exists(savePath)):
        os.makedirs(savePath)

    x = Positions[0]
    y = Positions[1]
    z = Positions[2]

    vx = Velocities[0]
    vy = Velocities[1]
    vz = Velocities[2]
    
    if (filetype == "txt"):

        #Header info: 
        info = "Plummer Sphere Model Random Initial Conditions from ic.py"
        info += "\n\n        masses             position_x           position_y           position_z           velocity_x           velocity_y           velocity_z"

        np.savetxt(savePath + filename + '.txt', list(zip(Masses, x, y, z, vx, vy, vz)), header = info, fmt="%20.8e")

    else:
        from snapwrite import write_snapshot

        pos_vectors = np.column_stack((x,y,z))
        vel_vectors = np.column_stack((vx,vy,vz))

        pos_vectors = np.array(pos_vectors, order='C')
        vel_vectors = np.array(vel_vectors, order='C')

        pos_vectors.shape = (1, -1)
        vel_vectors.shape = (1, -1)

        #Number of particles of each gadget type: 
        n_part = [0, 0, 0, 0, N, 0]   # only stars 

        #Gadget type of each particle
        IDs = np.arange(1, N + 1)

        data_list = [pos_vectors[0], vel_vectors[0], IDs, Masses]

        write_snapshot(n_part, data_list, outfile = savePath + filename, file_format = filetype)
        

if (__name__ == '__main__'):
    main()
