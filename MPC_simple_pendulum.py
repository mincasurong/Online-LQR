import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

class PendulumMPC:
    def __init__(self, params, Q, R, Np, Nc, dt, x0, u_min, u_max):
        self.g = params['g']
        self.l = params['l']
        self.m = params['m']
        self.Q = Q
        self.R = R
        self.Np = Np
        self.Nc = Nc
        self.dt = dt
        self.x = x0
        self.u_min = u_min
        self.u_max = u_max

    def dynamics(self, x, u):
        theta, theta_dot = x
        theta_ddot = -self.g / self.l * np.sin(theta) + u / (self.m * self.l**2)
        theta_next = theta + self.dt * theta_dot
        theta_dot_next = theta_dot + self.dt * theta_ddot
        return np.array([theta_next, theta_dot_next])

    def predict(self, u_seq):
        x_pred = self.x.copy()
        x_seq = [x_pred]
        for u in u_seq:
            x_pred = self.dynamics(x_pred, u)
            x_seq.append(x_pred)
        return np.array(x_seq)

    def cost_function(self, u_seq, x_ref):
        u_seq = np.reshape(u_seq, (self.Nc,))
        x_seq = self.predict(u_seq)
        cost = 0
        for k in range(min(self.Np, len(x_seq))):
            x_error = x_seq[k] - x_ref[k]
            if k < self.Nc:
                u = u_seq[k]
            else:
                u = 0
            cost += x_error.T @ self.Q @ x_error + u * self.R * u
        return cost

    def solve_mpc(self, x_ref):
        u_init = np.zeros((self.Nc,))
        bounds = [(self.u_min, self.u_max)] * self.Nc
        result = minimize(self.cost_function, u_init, args=(x_ref,), bounds=bounds)
        u_opt = result.x[0]
        self.x = self.dynamics(self.x, u_opt)
        return u_opt

# Pendulum parameters
params = {
    'g': 9.81,
    'l': 1.0,
    'm': 1.0
}

# MPC parameters
Q = np.diag([10, 1])
R = 0.1
Np = 20
Nc = 5
dt = 0.1
x0 = np.array([0, 0])  # Initial state: pendulum hanging down
u_min = -2
u_max = 2

# Reference trajectory: swing-up to the upright position
x_ref = np.zeros((Np + 1, 2))
x_ref[:, 0] = np.pi  # target angle

# Create MPC controller
mpc = PendulumMPC(params, Q, R, Np, Nc, dt, x0, u_min, u_max)

# Simulation
time = np.arange(0, 100, dt)
x_trajectory = [x0]
u_trajectory = []

for t in time:
    u = mpc.solve_mpc(x_ref)
    u_trajectory.append(u)
    x_trajectory.append(mpc.x)

x_trajectory = np.array(x_trajectory)
u_trajectory = np.array(u_trajectory)

# Plot results
plt.figure(figsize=(12, 6))
plt.subplot(3, 1, 1)
plt.plot(time, x_trajectory[:-1, 0])
plt.ylabel('Theta (rad)')
plt.title('Pendulum Swing-Up using MPC')

plt.subplot(3, 1, 2)
plt.plot(time, x_trajectory[:-1, 1])
plt.ylabel('Theta_dot (rad/s)')

plt.subplot(3, 1, 3)
plt.plot(time, u_trajectory)
plt.ylabel('Control Input (u)')
plt.xlabel('Time (s)')

plt.tight_layout()
plt.show()
