import matplotlib.pyplot as plt 
import numpy as np
from typing import List


def iterate(space, F):
    for i in range(1, space.shape[0]):
        for j in range(1, space.shape[1] - 1):
            space[i,j] = (1-2*F) * space[i-1,j] + F * space[i-1,j-1] + F * space[i-1,j+1]
            if (i + 1 == space.shape[0] and j + 1 == space.shape[1]):
                print(f"{(1-2*F) * space[i-1,j]} + {F * space[i-1,j-1]} + {F * space[i-1,j+1]}")
                print(f"{i}, {j}: {space[i,j]}")
                print(space[i-10:i+1, :])
                print("============")

def C_n(n : int):
    if n == 1:
        return -4 / np.pi
    else:
        return (1 - np.cos((n+1)*np.pi))/ (n+1) / np.pi - (1-np.cos((1-n)*np.pi)) / (1-n)/np.pi + 4*np.cos(np.pi*n) / n / np.pi

def f_x(x : np.array, t : np.array, order : int) -> float:
    if order < 1:
        raise Exception("Order must be at least 2")
    # Shape: (len(t), len(x))
    values = 2 * np.tile(x, (t.shape[0], 1))
    for n in range(1, order + 1):
        coeff = C_n(n)
        sin_term = np.sin(np.pi * n * x) # shape: (len(x),)
        decay = np.exp(-2 * (np.pi * n)**2 * t[:, None]) # shape: (len(t), 1)
        values += coeff * sin_term[None, :] * decay # broadcasted multiplication
    return values # shape: (len(t), len(x))

def analytical_solver(times: List[float], order: int = 100, num_points: int = 100):
    pass  # Keep this function if needed for other purposes


times = np.array([0, 0.1, 0.2])  # Example values for times
num_points = 100  # Define num_points
x_points = np.linspace(0, 1, num_points)

order = 100  # Define order before using it
values = f_x(x_points, times, order)

plt.figure(figsize=(10, 6))
for i, t in enumerate(times):
    plt.plot(x_points, values[i, :], label=f"t = {t}")
    plt.plot(x_points, values[i, :], label=f"t = {t}")

plt.xlabel("x")
plt.ylabel("u(x,t)")
plt.title(f"Analytical Solution of the Heat Equation(N={order})")  # Ensure order is defined
plt.title(f"Analytical Solution of the Heat Equation(N={order})")
plt.legend(loc="upper left")
plt.grid(True)
plt.savefig("plots/2a_analytical.png", dpi=300)
plt.show()