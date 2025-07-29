import numpy as np
import matplotlib.pyplot as plt

def getConserved(rho, vx, vy, P, gamma, vol):
    Mass = rho * vol
    Momx = Mass * vx
    Momy = Mass * vy
    Energy = P / (gamma - 1) * vol + 0.5 * (Momx**2 + Momy**2) / Mass
    return Mass, Momx, Momy, Energy


def getPrimitive(Mass,Momx,Momy,Energy,gamma,vol):
    rho = Mass / vol
    vx = Momx / Mass
    vy = Momy / Mass
    P = (Energy / vol - 0.5 * rho * (vx**2 + vy**2)) *(gamma - 1)   
    return rho,vx,vy,P


def getGradient(f, dx):
    R = -1
    L = 1
    f_dx = ( np.roll(f, R, axis=0) - np.roll(f, L, axis=0) ) / (2 * dx)
    f_dy = ( np.roll(f, R, axis=1) - np.roll(f, L, axis=1) ) / (2 * dx)
    return f_dx, f_dy

def main():

    #control param
    N           = 128 # the number of cells
    box_size   = 1.0 # the size of the box
    gamma       = 1.4 # the ratio of specific heats
    t_output    = 0.1 # the time interval of output
    useSlopeLimiter = False # whether to use slope limiter
    plotRealTime = False # whether to plot in real time
    t           = 0   # the current time
    t_end       = 2 # the end time of simulation
    CFL_number  = 0.5 # the coutant factor number

    #Mesh param
    dx          = box_size / N
    vol         = dx**2*1
    x_c_list    = np.linspace(0.5*dx,box_size-0.5*dx,N)# the center of the cell
    Y,X         = np.meshgrid(x_c_list,x_c_list)#  x(y)crood at [i,j]
    
    print("===== ====== ====== ====== ===== ===== ===== ===== ===== ===== =====")
    print(X)
    print("===== ====== ====== ====== ===== ===== ===== ===== ===== ===== =====")
    print(Y)
    print("===== ====== ====== ====== ===== ===== ===== ===== ===== ===== =====")

    #initial Conditions
    w0 = 0.1
    sigma = 0.05 / np.sqrt(0.5*np.pi)
    rho = 1 + w0 * np.exp(-((X - 0.5)**2 + (Y - 0.5)**2) / (2 * sigma**2))
    vx  = np.zeros_like(rho)
    vy  = np.zeros_like(rho)
    P = (gamma - 1) * (0.5 * (rho * (vx**2 + vy**2)))

    #Get the conserved variables
    Mass, Momx, Momy = getConserved(rho,vx,vy,P,gamma,vol)
    
    #plot figure
    fig = plt.figure(figsize =(4,4),dpi = 100)
    outputCount = 1

    #simulation loop
    while t < t_end:
        # get the Primitive var
        rho,vx,vy,p = getPrimitive(Mass,Momx,Momy,Energy,gamma,vol)
        
        #get the time step, dt = CFL * min(dx / (a + v)), where a is the sound speed
        dt = CFL_number * np.min(dx / (np.sqrt(gamma * p / rho)) + np.sqrt(vx**2 + vy**2))

        plotThisTurn= False# whether to plot this turn
        if t + dt > t_output * outputCount:
            dt = t_output * outputCount - t
            plotThisTurn = True
        
        # cal the gradient
        rho_dx,rho_dy = getGradient(rho, dx)
        vx_dx,vx_dy = getGradient(vx, dx)
        vy_dx,vy_dy = getGradient(vy, dx)
        P_dx,P_dy = getGradient(P, dx)
        
        # slope limiter
        if useSlopeLimiter:
            rho_dx,rho_dy = slopeLimiter(rho,dx,rho_dx,rho_dy)
            vx_dx,vx_dy = slopeLimiter(vx,dx,vx_dx,vx_dy)
            vy_dx,vy_dy = slopeLimiter(vy,dx,vy_dx,vy_dy)
            P_dx,P_dy = slopeLimiter(P,dx,P_dx,P_dy)
        
        # use the MUSCL scheme to calculate the flux
        rho_prime = rho - 0.5 * dt * (rho_dx * vx + rho_dy * vy + rho * (vx_dx + vy_dy))
        vx_prime = vx - 0.5 * dt * (vx_dx * vx + vy_dx * vy + P_dx / rho)
        vy_prime = vy - 0.5 * dt * (vx_dy * vx + vy_dy * vy + P_dy / rho)
        P_prime = P - 0.5 * dt * (P_dx * vx + P_dy * vy + gamma * P / rho * (vx_dx + vy_dy))
        
        #extrapolate in space
        rho_XL,rho_XR, rho_YL, rho_YR = extrapolate(rho_prime,rho_dx,rho_dy, dx)
        vx_XL,vx_XR, vx_YL, vx_YR = extrapolate(vx_prime,vx_dx,vx_dy, dx)
        vy_XL,vy_XR, vy_YL, vy_YR = extrapolate(vy_prime,vy_dx,vy_dy, dx)
        P_XL,P_XR, P_YL, P_YR = extrapolate(P_prime,P_dx,P_dy, dx)
        
        # calculate the flux
        flux_mass_x, flux_momx_x,flux_momy_x,flux_energy_x = getFlux(rho_XL,rho_XR,vx_XL,vx_XR,vy_XL,vy_XR,P_XL,P_XR,gamma)
        flux_mass_y, flux_momx_y,flux_momy_y,flux_energy_y = getFlux(rho_YL,rho_YR,vx_YL,vx_YR,vy_YL,vy_YR,P_YL,P_YR,gamma)
        
        # update the solution
        Mass = applyFlux(Mass, flux_mass_x, flux_mass_y, dx, dt)
        Momx = applyFlux(Momx, flux_momx_x, flux_momx_y, dx, dt)
        Momy = applyFlux(Momy, flux_momy_x, flux_momy_y, dx, dt)
        Energy = applyFlux(Energy, flux_energy_x, flux_energy_y, dx, dt)
        
        # update the time
        t += dt
        
        #plot the result
        if(plotRealTime and plotThisTurn) or (t>=t_end):
            plt.cla()
            plt.imshow(rho.T)
            plt.clim(0.8,2.2)
            ax = plt.gca()
            ax.invert_yaxis()
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.set_aspect('equal', adjustable='box')
            plt.pause(0.001)
            outputCount += 1
            
    # save the final result
    plt.savefig("FVM-2D.png", dpi=300)
    plt.show()
if __name__=="__main__":
    main()
