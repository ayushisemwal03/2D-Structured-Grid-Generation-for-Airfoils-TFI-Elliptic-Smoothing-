#   IMPORTING REQUIRED MODULES
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

def NACA_AIRFOIL(file_path):
    raw_data = pd.read_csv(file_path, skiprows=5, header=None, names=["x", "y"])

    # Remove non-numeric rows and convert to float
    airfoil_coords = raw_data[pd.to_numeric(raw_data["x"], errors="coerce").notnull()]
    airfoil_coords["x"] = airfoil_coords["x"].astype(float)/1000
    airfoil_coords["y"] = airfoil_coords["y"].astype(float)/1000

    if not np.allclose([airfoil_coords.iloc[0]["x"], airfoil_coords.iloc[0]["y"]],
                       [airfoil_coords.iloc[-1]["x"], airfoil_coords.iloc[-1]["y"]]):
        airfoil_coords = pd.concat([airfoil_coords, airfoil_coords.iloc[[0]]], ignore_index=True)

    return airfoil_coords["x"].values, airfoil_coords["y"].values

def CUBICSPLINE_FITTING(x, y, n_points):    #   Spline Fitting for smoothness
    t = np.linspace(0, 1, len(x))
    cs_x, cs_y = CubicSpline(t, x), CubicSpline(t, y)
    t_new = np.linspace(0, 1, n_points)
    return cs_x(t_new), cs_y(t_new)


#   Function to have grid using TFI

def TFI_GRID_GENERATION(NX, NY, chord_length, radius):
    filepath = "NACA63412coordinates.csv"
    x_airfoil, y_airfoil = NACA_AIRFOIL(filepath)
    x_airfoil, y_airfoil = CUBICSPLINE_FITTING(x_airfoil, y_airfoil, NX)
    
    delta_xi, delta_eta = 1 / (NX - 1), 1 / (NY - 1)    #   Step Length
    xi_val, eta_values = np.linspace(0, 1, NX), np.linspace(0, 1, NY)

    #   Corner Points of Computational Domai

    X_A,Y_A = chord_length / 2 , 0
    X_B,Y_B = chord_length / 2 , 0
    X_C,Y_C = radius , 0
    X_D,Y_D = radius , 0
    
    #   Boundary Conditions

    x_bottom, y_bottom = x_airfoil, y_airfoil

    x_top , y_top = radius * np.cos(2 * np.pi * xi_val), radius * np.sin(2 * np.pi * xi_val)

    x_left , y_left = (1 - eta_values) * X_A + eta_values * X_D , (1 - eta_values) * Y_A + eta_values * Y_D
    x_right , y_right = (1 - eta_values) * X_B + eta_values * X_C , (1 - eta_values) * Y_B + eta_values * Y_C
    
    X, Y = np.zeros((NX, NY)), np.zeros((NX, NY))
    
    for i in range(NX):             #   The Formula can be Optimized does not go for it, to have same look of formula as in Notes.
        xi = i * delta_xi
        for j in range(NY):
            eta = j * delta_eta
            X[i, j] = (1 - eta) * x_bottom[i] + eta * x_top[i] + (1 - xi) * x_left[j] + xi * x_right[j] - (
                (1 - xi) * (1 - eta) * X_A + xi * (1 - eta) * X_B + (1 - xi) * eta * X_D + xi * eta * X_C)
            
            Y[i, j] = (1 - eta) * y_bottom[i] + eta * y_top[i] + (1 - xi) * y_left[j] + xi * y_right[j] - (
                (1 - xi) * (1 - eta) * Y_A + xi * (1 - eta) * Y_B + (1 - xi) * eta * Y_D + xi * eta * Y_C)
    
    return X, Y, x_airfoil, y_airfoil

#   Coefficients for Elliptic Grid Generation Method

def abc(X, Y):
    x_xi, x_eta = np.gradient(X, axis=0), np.gradient(X, axis=1)
    y_xi, y_eta = np.gradient(Y, axis=0), np.gradient(Y, axis=1)
    return x_eta**2 + y_eta**2, x_xi * x_eta + y_xi * y_eta, x_xi**2 + y_xi**2

#   Function to have grid using ELLIPTIC GRID GENERATION METHOD

def ELLIPTIC_GRID_GENERATION(X, Y, iterations=5000, tol=1e-3):
    X_new, Y_new = X.copy(), Y.copy()
    delta_xi, delta_eta = 1.0 / (NX - 1), 1.0 / (NY - 1)

    for _ in range(iterations):
        X_old, Y_old = X_new.copy(), Y_new.copy()
        a, b, c = abc(X_old, Y_old)

        for i in range(1, NX-1):
            for j in range(1, NY-1):
                X_new[i, j] = (a[i, j] * (X_old[i+1, j] + X_old[i-1, j]) / delta_xi**2 +c[i, j] * (X_old[i, j+1] + X_old[i, j-1]) / delta_eta**2 -b[i, j] * (X_old[i+1, j+1] - X_old[i+1, j-1] - X_old[i-1, j+1] + X_old[i-1, j-1]) / (2 * delta_xi * delta_eta)) /(2 * (a[i, j] / delta_xi**2 + c[i, j] / delta_eta**2))
                
                Y_new[i, j] = (a[i, j] * (Y_old[i+1, j] + Y_old[i-1, j]) / delta_xi**2 +c[i, j] * (Y_old[i, j+1] + Y_old[i, j-1]) / delta_eta**2 -b[i, j] * (Y_old[i+1, j+1] - Y_old[i+1, j-1] - Y_old[i-1, j+1] + Y_old[i-1, j-1]) / (2 * delta_xi * delta_eta)) /(2 * (a[i, j] / delta_xi**2 + c[i, j] / delta_eta**2))
        
        for j in range(1, NY-1):
                X_new[0, j] = (a[0, j] * (X_old[1, j] + X_old[-2, j]) / delta_xi**2 +c[0, j] * (X_old[0, j+1] + X_old[0, j-1]) / delta_eta**2 -b[0, j] * (X_old[1, j+1] - X_old[1, j-1] - X_old[-2, j+1] + X_old[-2, j-1]) / (2 * delta_xi * delta_eta)) /(2 * (a[0, j] / delta_xi**2 + c[0, j] / delta_eta**2))
                
                Y_new[0, j] = (a[0, j] * (Y_old[1, j] + Y_old[-2, j]) / delta_xi**2 +c[0, j] * (Y_old[0, j+1] + Y_old[0, j-1]) / delta_eta**2 -b[0, j] * (Y_old[1, j+1] - Y_old[1, j-1] - Y_old[-2, j+1] + Y_old[-2, j-1]) / (2 * delta_xi * delta_eta)) /(2 * (a[0, j] / delta_xi**2 + c[0, j] / delta_eta**2))
        
        X_new[-1,:] = X_new[0,:]
        Y_new[-1,:] = Y_new[0,:]
       
        error = np.max(np.abs(X_new - X_old)) + np.max(np.abs(Y_new - Y_old))
        print(f'Iteration {_}, Error: {error}')
        if error < tol:
            break
        
    return X_new, Y_new

def plot_grid(X, Y, x_airfoil, y_airfoil, title):
    plt.figure(figsize=(10, 8))
    plt.plot(X, Y, linewidth=1)
    plt.plot(X.T, Y.T, linewidth=1)
    plt.plot(x_airfoil, y_airfoil, color='red', linewidth=2)
    plt.xlabel("x"),plt.ylabel("y"),plt.title(title),plt.axis("equal")
    plt.show()

# PLOTTING

NX, NY = 60, 60
X_TFI, Y_TFI, x_airfoil, y_airfoil = TFI_GRID_GENERATION(NX, NY, 1, 10)
plot_grid(X_TFI, Y_TFI, x_airfoil, y_airfoil, "Grid Generation using TFI")

X_ELLIPTIC, Y_ELLIPTIC = ELLIPTIC_GRID_GENERATION(X_TFI, Y_TFI)
plot_grid(X_ELLIPTIC, Y_ELLIPTIC, x_airfoil, y_airfoil, "Grid Generation using Elliptic Solver")

