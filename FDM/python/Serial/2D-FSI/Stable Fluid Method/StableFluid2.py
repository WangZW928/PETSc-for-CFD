import numpy as np
import scipy.sparse.linalg as splinalg
from scipy import interpolate
import matplotlib.pyplot as plt
import math

# Optional
import cmasher as cmr
from tqdm import tqdm


def forcing_function(time, point):
    
    
    if(time>10):
        #F = 1.0
        F = 0
    else:
        F = 0
    theta = math.sin((time-10)/5) * math.pi / 6
    fx = math.cos(theta) * F
    fy = math.sin(theta) * F

    return forced_value


def main():

    forcing_function_vectorized = np.vectorize(
        pyfunc=forcing_function,
        signature="(),(d)->(d)",
    )

    

    


    
    
    
    
    
    def curl_2d(vector_field):
        curl_applied = (
            partial_derivative_x(vector_field[..., 1])
            -
            partial_derivative_y(vector_field[..., 0])
        )

        return curl_applied

    def advect(field, vector_field):
        backtraced_positions = np.clip(
            (
                coordinates
                -
                TIME_STEP_LENGTH
                *
                vector_field
            ),
            0.0,
            DOMAIN_SIZE,
        )

        advected_field = interpolate.interpn(
            points=(x, y),
            values=field,
            xi=backtraced_positions,
        )

        return advected_field
    
# 半拉格朗日对流
# u(t_n,x) = u(t_(n-1),x_l)
# x_L = x - timestep * u(t_(n-1),x_l)
    
    def diffusion_operator(vector_field_flattened):
        vector_field = vector_field_flattened.reshape(vector_shape)

        diffusion_applied = (
            vector_field
            -
            KINEMATIC_VISCOSITY
            *
            TIME_STEP_LENGTH
            *
            laplace(vector_field)
        )

        return diffusion_applied.flatten()
    
    def poisson_operator(field_flattened):
        field = field_flattened.reshape(scalar_shape)

        poisson_applied = laplace(field)

        return poisson_applied.flatten()

    plt.style.use("dark_background")
    plt.figure(figsize=(5, 5), dpi=160)

    velocities_prev = np.zeros(vector_shape)
    #velocities_prev[:, :, 0] = 1.0  # u分量
    
    time_current = 0.0
    for i in tqdm(range(N_TIME_STEPS)):
        time_current += TIME_STEP_LENGTH
        
        velocities_prev[0, :, 0] = 1.0  # u分量
        velocities_prev[0, :, 1] = 0.0  # v分量
        velocities_prev[-1, :, :] = velocities_prev[-2, :, :]

        velocities_prev[:, -1, :] =  velocities_prev[:, -2, :]

        velocities_prev[:, 0, :] = velocities_prev[:, 1, :]

        

        forces = forcing_function_vectorized(
            time_current,
            coordinates,
        )

        # (1) Apply Forces
        velocities_forces_applied = (
            velocities_prev
            +
            TIME_STEP_LENGTH
            *
            forces
        )

        # (2) Nonlinear convection (=self-advection)
        velocities_advected = advect(
            field=velocities_forces_applied,
            vector_field=velocities_forces_applied,
        )

        # (3) Diffuse
        velocities_diffused = splinalg.cg(
            A=splinalg.LinearOperator(
                shape=(vector_dof, vector_dof),
                matvec=diffusion_operator,
            ),
            b=velocities_advected.flatten(),
            maxiter=MAX_ITER_CG,
        )[0].reshape(vector_shape)

        # (4.1) Compute a pressure correction
        pressure = splinalg.cg(
            A=splinalg.LinearOperator(
                shape=(scalar_dof, scalar_dof),
                matvec=poisson_operator,
            ),
            b=divergence(velocities_diffused).flatten(),
            maxiter=MAX_ITER_CG,
        )[0].reshape(scalar_shape)

        # (4.2) Correct the velocities to be incompressible
        velocities_projected = (
            velocities_diffused
            -
            gradient(pressure)
        )

        # Advance to next time step
        velocities_prev = velocities_projected

        # Plot
        curl = curl_2d(velocities_projected)
        plt.contourf(
            X,
            Y,
            curl,
            cmap=cmr.redshift,
            levels=100,
        )
        plt.quiver(
            X, 
            Y,
            velocities_projected[..., 0],
            velocities_projected[..., 1],
            color="dimgray",
        )
        plt.draw()
        plt.pause(0.0001)
        plt.clf()

    plt.show()
        


if __name__ == "__main__":
    main()