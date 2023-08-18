import fym
import matplotlib.pyplot as plt
import numpy as np
from fym.utils.rot import angle2quat
from scipy.integrate import quad
from scipy.optimize import minimize

from simpledyn import simpleLC62

class TrajOpt(fym.BaseEnv):
    ENV_CONFIG = {
        "fkw": {
            "dt": 0.01,
            "max_t": 10,
        },
        "plant": {
            "init": {
                "pos": np.vstack((0, 0, 0)),
                "vel": np.vstack((0, 0, 0)),
            },
        },
    }

    

    def __init__(self, env_config={}):
        env_config = safeupdate(self.ENV_CONFIG, env_config)
        super().__init__(**env_config["fkw"])
        self.plant = simpleLC62(env_config["plant"])
        self.tf = env_config["fkw"]["max_t"]
        self.cruise_height = 10
        self.cruise_speed = 45


    def step(self):
        env_info, done = self.update()
        return done, env_info

    def observation(self):
        return self.observe_flat()

    
    def set_dot(self, t):

        # TODO
        # U = 
        self.plant.set_dot(t, U)


        env_info = {
            "t": t,
            **self.observe_dict(),
            "U": U,
        }

        return env_info

    def minimize_cost(
        self,
        tf,
        U0={"F_r": 0, "F_p": 0, "theta": 0},
        method="SLSQP",
        options={"disp": False, "ftol": 1e-10},
    ):

        U0 = list(U0.values())
        result = minimize(
            self.performance_index,
            U0,
            args=(tf,),
            method=method,
            options=options,
        )
        
        U = np.vstack((result.x))
        return U

    def performance_index(self, U, tf):
        energy = self.energy(U)
        J = quad(energy, 0, tf)
        return J

    def energy(self, U):

       pos, vel = self.observe_list()
       Vx, Vz = np.ravel(vel)

       F_r, F_p, theta = np.ravel(U)

       P_r = F_r * Vz
       P_p = F_p * Vx
       P = P_r + P_p
       return P
    
   def Boundary_cond(
       self,
       fixed={"h", "V"},
   ):

       # TODO
       # 
       return X_trim, U_trim
