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

        tf = self.tf
        U = self.minimize_cost(tf)
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
    

def run():
    env = TrajOpt()
    flogger = fym.Logger("data.h5")

    env.reset()
    try:
        while True:
            env.render()

            done, env_info = env.step()
            flogger.record(env=env_info)

            if done:
                break

    finally:
        flogger.close()
        plot()


def plot():
    data = fym.load("data.h5")["env"]

    """ Figure 1 - States """
    fig, axes = plt.subplots(2, 2, squeeze=False, sharex=True)

    """ Column 1 - States: Position """
    ax = axes[0, 0]
    ax.plot(data["t"], data["plant"]["pos"][:, 0].squeeze(-1), "k-")
    ax.set_ylabel(r"$x$, m")
    ax.set_xlabel("Time, sec")
    ax.set_xlim(data["t"][0], data["t"][-1])

    ax = axes[1, 0]
    ax.plot(data["t"], data["plant"]["pos"][:, 1].squeeze(-1), "k-")
    ax.set_ylabel(r"$y$, m")
    ax.set_xlabel("Time, sec")
    ax.set_xlim(data["t"][0], data["t"][-1])

    """ Column 2 - States: Velocity """
    ax = axes[0, 1]
    ax.plot(data["t"], data["plant"]["vel"][:, 0].squeeze(-1), "k-")
    ax.set_ylabel(r"$Vx$, m")
    ax.set_xlabel("Time, sec")
    ax.set_xlim(data["t"][0], data["t"][-1])

    ax = axes[1, 1]
    ax.plot(data["t"], data["plant"]["vel"][:, 1].squeeze(-1), "k-")
    ax.set_ylabel(r"$Vy$, m")
    ax.set_xlabel("Time, sec")
    ax.set_xlim(data["t"][0], data["t"][-1])


    """ Figure 2 - Optimization Variables """
    fig, axes = plt.subplots(3, 1, squeeze=False, sharex=True)

    ax = axes[0, 0]
    ax.plot(data["t"], data["U"][:, 0].squeeze(-1), "k-")
    ax.set_ylabel(r"$F_r$, N")
    ax.set_xlabel("Time, sec")
    ax.set_xlim(data["t"][0], data["t"][-1])

    ax = axes[1, 0]
    ax.plot(data["t"], data["U"][:, 1].squeeze(-1), "k-")
    ax.set_ylabel(r"$F_p$, N")
    ax.set_xlabel("Time, sec")
    ax.set_xlim(data["t"][0], data["t"][-1])

    ax = axes[2, 0]
    ax.plot(data["t"], data["U"][:, 2].squeeze(-1), "k-")
    ax.set_ylabel(r"$\theta$, deg")
    ax.set_xlabel("Time, sec")
    ax.set_xlim(data["t"][0], data["t"][-1])

    
    plt.show()


def main(args):
    if args.only_plot:
        plot()
        return
    else:
        run()
        if args.plot:
            plot()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--plot", action="store_true")
    parser.add_argument("-P", "--only-plot", action="store_true")
    args = parser.parse_args()
    main(args)

