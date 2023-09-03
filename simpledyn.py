"""
Simple Longitudinal Dynamic Model of LC62-50B
"""

import fym
import numpy as np
import scipy
from numpy import cos, sin


class simpleLC62(fym.BaseEnv):

    g = 9.81
    m = 41.97

    S1 = 0.2624  # [m^2]
    S2 = 0.5898  # [m^2]
    S = S1 + S2
    St = 0.06894022  # [m^2]
    c = 0.551  # Main wing chord length [m]
    b = 1.1  # Main wing half span [m]
    d = 0.849  # Moment arm length [m]

    tables = {
        "alp": np.deg2rad(np.array([0, 2, 4, 6, 8, 10, 12, 14, 16])),  # [rad]
        "CL": np.array(
            [0.1931, 0.4075, 0.6112, 0.7939, 0.9270, 1.0775, 0.9577, 1.0497, 1.0635]
        ),
        "CD": np.array(
            [0.0617, 0.0668, 0.0788, 0.0948, 0.1199, 0.1504, 0.2105, 0.2594, 0.3128]
        ),
        "Cm": np.array(
            [
                0.0406,
                0.0141,
                -0.0208,
                -0.0480,
                -0.2717,
                -0.4096,
                -0.1448,
                -0.2067,
                -0.2548
            ]
        ),
        "cmd": np.array(
            [
                0,
                0.2,
                0.255,
                0.310,
                0.365,
                0.420,
                0.475,
                0.530,
                0.585,
                0.640,
                0.695,
                0.750,
            ]
        ),
        "th_p": np.array(
            [
                0,
                1.39,
                4.22,
                7.89,
                12.36,
                17.60,
                23.19,
                29.99,
                39.09,
                46.14,
                52.67,
                59.69,
            ]
        ),
        "tq_p": np.array(
            [
                0, 
                0.12, 
                0.35, 
                0.66, 
                1.04, 
                1.47, 
                1.93, 
                2.50, 
                3.25, 
                3.83, 
                4.35, 
                4.95,
            ]
        ),
        "th_r": [-19281, 36503, -992.75, 0],
        "tq_r": [-6.3961, 12.092, -0.3156, 0],
    }

    control_limits = {
        "cmd": (0, 1),
        "dela": np.deg2rad((-10, 10)),
        "dele": np.deg2rad((-10, 10)),
        "delr": np.deg2rad((-10, 10)),
    }

    ENV_CONFIG = {
        "init": {
            "pos": np.zeros((2, 1)),
            "vel": np.zeros((2, 1)),
        },
    }

    def __init__(self, env_cinfig={}):
        env_config = safeupdate(self.ENV_CONFIG, env_config)
        super().__init__()
        self.pos = fym.BaseSystem(env_config["init"]["pos"])
        self.vel = fym.BaseSystem(env_config["init"]["vel"])


    def set_dot(self, t, U):
        pos, vel = self.observe_list()
        F_w = self.get_Fw(pos[1], vel)
        Fx_w, Fz_w = np.ravel(F_w)

        # TODO 
        # find U(t)
        F_r, F_p, theta = np.ravel(U)

        R = np.array(
            [
                [cos(theta), -sin(theta)], 
                [sin(theta),  cos(theta)]
            ]
        )

        Vx_dot = (F_p + Fx_w)/self.m - self.g * sin(theta)
        Vz_dot = (-F_r + Fz_w)/self.m + self.g * cos(theta)
      
        self.pos.dot = R @ vel
        self.vel.dot = np.vstack((Vx_dot, Vz_dot))


    def get_Fw(self, z, vel):
        Vx, Vz = np.ravel(vel)
        V = np.linalg.norm(vel)

        rho = self.get_rho(-z) 
        qbar = 0.5 * rho * V**2
        alp = np.arctan2(Vz, Vx)
        CL, CD, _ = self.aero_coeff(alp)
        S = self.S

        Fx_w = -qbar * S * CD
        Fz_w = -qbar * S * CL

        F_w = np.vstack((Fx_w, Fz_w))
        return F_w


    def get_rho(self, altitude):
        pressure = 101325 * (1 - 2.25569e-5 * altitude) ** 5.25616
        temperature = 288.14 - 0.00649 * altitude
        return pressure / (287 * temperature)


    def aero_coeff(self, alp):
        CL = np.interp(alp, self.tables["alp"], self.tables["CL"])
        CD = np.interp(alp, self.tables["alp"], self.tables["CD"])
        Cm = np.interp(alp, self.tables["alp"], self.tables["Cm"])
        return np.vstack((CL, CD, Cm))


    def get_trim(
        self,
        q0={
            "alp": 0.0,
            "F_r": 0,
            "F_p": 0,
            "theta": 0,
        },
        fixed={"h": 10, "V": 45},
        method="SLSQP",
        options={"disp": False, "ftol": 1e-10},
    ):
        q0 = list(q0.values())
        fixed = list(fixed.values())

        Fr_max = np.polyval(self.tables["th_r"], 1) * self.g / 1000
        Ft_max = np.interp(1, self.tables["cmd"], self.tables["th_p"])
        
        bounds = (
            np.deg2rad((0, 20)), 
            ((0, Fr_max)),
            ((0, Ft_max)),
            ((0, np.pi / 2)),
        )

        result = scipy.optimize.minimize(
            self._trim_cost,
            q0,
            args=(fixed,),
            bounds=bounds,
            method=method,
            options=options,
        )

        alp_trim, Fr_trim, Ft_trim, theta_trim = result.x

        z_trim, V_trim = fixed
        Vx_trim = V_trim * cos(alp_trim)
        Vz_trim = V_trim * sin(alp_trim)
        vel_trim = np.vstack((Vx_trim, Vz_trim))

        X_trim = (z_trim, vel_trim)
        U_trim = (Fr_trim, Ft_trim, theta_trim)

        return X_trim, U_trim

    def _trim_cost(self, q, fixed):
        h, V = fixed
        alp, F_r, F_p, theta = q

        z = h
        vel = V * np.vstack((cos(alp), sin(alp)))
        F_w = self.get_Fw(h, vel)
        Fx_w, Fz_w = np.ravel(F_w)
        Vx_dot = (F_p + Fx_w) / self.m - self.g * sin(theta)
        Vz_dot = (-F_r + Fz_w) / self.m + self.g * cos(theta)
        vel_dot = np.vstack((Vx_dot, Vz_dot))

        return vel_dot.T * vel_dot


