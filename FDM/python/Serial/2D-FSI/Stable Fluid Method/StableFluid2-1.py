import numpy as np
import scipy.sparse.linalg as splinalg
from scipy import interpolate
import matplotlib.pyplot as plt
import math

# Optional
import cmasher as cmr
from tqdm import tqdm

class StableFluid2:
    DOMAIN_SIZE = 1.0
    N_POINTS = 41
    N_TIME_STEPS = 2000
    TIME_STEP_LENGTH = 0.1
    KINEMATIC_VISCOSITY = 0.0001
    MAX_ITER_CG = None

    def __init__(self):
        
        self.lx = self.DOMAIN_SIZE
        self.ly = self.DOMAIN_SIZE
        self.nx = self.N_POINTS
        self.dx = self.lx / (self.nx - 1)
        self.ny = self.N_POINTS
        self.nt = self.N_TIME_STEPS
        self.dt = self.TIME_STEP_LENGTH 
        self.time = 0
        self.nu = self.KINEMATIC_VISCOSITY
        self.max_iter_cg = self.MAX_ITER_CG
        
        self.scalar_shape = (self.nx, self.nx)
        self.scalar_dof = self.nx * self.ny
        self.vector_shape = (self.nx, self.ny, 2)
        self.vector_dof = 2 * self.scalar_dof
        
        self.x = np.linspace(0, self.lx, self.nx)
        self.y = np.linspace(0, self.ly, self.ny)
        self.X, self.Y = np.meshgrid(self.x, self.y,indexing="ij")
        
        self.coordinates = np.concatenate(
        (
            self.X[..., np.newaxis],
            self.Y[..., np.newaxis],
        ),    
        axis=-1,
        )
        
        self.Turbine_D = 0.2 * self.lx
        self.Turbine_R = 0.5 * self.Turbine_D
        self.Turbine_Center_x = 0.3 * self.lx
        self.Turbine_Center_y = 0.5 * self.ly
        self.Turbine_Center_x_index = int(self.Turbine_Center_x / self.dx)
        self.Turbine_Center_y_index = int(self.Turbine_Center_y / self.dy)
        self.Turbine_Inflow_x_index = int((self.Turbine_Center_D) / self.dx)
        self.Turbine_A = 0.25 * math.pi * (self.Turbine_R**2)
        self.Turbine_a = 0.25
        self.Turbine_Ct = 4 * self.Turbine_a * (1 - self.Turbine_a)
        self.Turbine_ud = 0
        self.rho = 1.0
    def forcing_function(self,ud,point):
        F = 0.5 * self.rho * self.A * self.Ct * (self.ud**2) / (1-self.Turbine_a)**2
        turbine_x = self.Turbine_Center_x_index * self.dx
        y_s = self.Turbine_Center_y - self.Turbine_R
        y_e = self.Turbine_Center_y + self.Turbine_R
        fx = 1
        fy = 0
        
        forced_value = (
        F
        *
        np.where(
            (
                (point[0] == turbine_x)
                &
                (point[1] > y_s)
                &
                (point[1] < y_e)
            ),
            np.array([-fx, fy]),
            np.array([0.0, 0.0]),
        )
    )

        return forced_value
    
    def get_ud(self):
        u = self.u
        self.Turbine_ud = 0  
        
    def partial_derivative_x(self):
        diff = np.zeros_like(self.field)

        diff[1:-1, 1:-1] = (
            (
                self.field[2:  , 1:-1]
                -
                self.field[0:-2, 1:-1]
            ) / (
                2 * self.dx
            )
        )

        return diff
    
    def partial_derivative_y(self):
        diff = np.zeros_like(self.field)

        diff[1:-1, 1:-1] = (
            (
                self.field[1:-1, 2:  ]
                -
                self.field[1:-1, 0:-2]
            ) / (
                2 * self.dy
            )
        )

        return diff  
    
    def laplace(self):
        diff = np.zeros_like(self.field)

        diff[1:-1, 1:-1] = (
            (
                self.field[0:-2, 1:-1]
                +
                self.field[1:-1, 0:-2]
                - 4 *
                self.field[1:-1, 1:-1]
                +
                self.field[2:  , 1:-1]
                +
                self.field[1:-1, 2:  ]
            ) / (
                self.dx**2
            )
        )

        return diff
    
    def divergence(self):
        divergence_applied = (
            self.partial_derivative_x(self.vector_field[..., 0])
            +
            self.partial_derivative_y(self.vector_field[..., 1])
        )

        return divergence_applied
    def gradient(self):
        gradient_applied = np.concatenate(
            (
                self.partial_derivative_x(self.field)[..., np.newaxis],
                self.partial_derivative_y(self.field)[..., np.newaxis],
            ),
            axis=-1,
        )

        return gradient_applied
