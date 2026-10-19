import numpy as np

# thanks to a Google search 2026/10/18 (AI generated)

def lambert_solver(r1, r2, tof, mu, lw=False, max_iter=100, tol=1e-12):
    """
    Solves Lambert's problem using Izzo's algorithm.
    
    Parameters:
        r1 (array_like): Position vector at departure (shape: 3,)
        r2 (array_like): Position vector at arrival (shape: 3,)
        tof (float): Time of flight
        mu (float): Gravitational parameter of the central body
        lw (bool): Long way flag. True for > 180 deg transfer, False for < 180 deg.
        max_iter (int): Maximum number of iterations for the Secant solver.
        tol (float): Convergence tolerance.
        
    Returns:
        v1 (np.ndarray): Velocity vector at departure (shape: 3,)
        v2 (np.ndarray): Velocity vector at arrival (shape: 3,)
    """
    r1 = np.asarray(r1, dtype=np.float64)
    r2 = np.asarray(r2, dtype=np.float64)
    
    R1 = np.linalg.norm(r1)
    R2 = np.linalg.norm(r2)
    
    # 1. Calculate chord and semi-perimeter
    c = np.linalg.norm(r2 - r1)
    s = 0.5 * (R1 + R2 + c)
    
    # 2. Determine transfer angle geometry
    cross_prod = np.cross(r1, r2)
    # Check if transfer is long-way or short-way
    if (cross_prod[2] >= 0 and not lw) or (cross_prod[2] < 0 and lw):
        ir1 = r1 / R1
        ir2 = r2 / R2
        ih = np.cross(ir1, ir2)
        ih /= np.linalg.norm(ih)
    else:
        ir1 = r1 / R1
        ir2 = r2 / R2
        ih = np.cross(ir2, ir1)
        ih /= np.linalg.norm(ih)
        
    it1 = np.cross(ih, ir1)
    it2 = np.cross(ih, ir2)
    
    # 3. Non-dimensionalize parameters
    lambda_param = np.sqrt(1.0 - c / s)
    if cross_prod[2] < 0:
        lambda_param = -lambda_param if lw else lambda_param
    else:
        lambda_param = -lambda_param if lw else lambda_param
        
    # Re-evaluate lambda sign based on short/long way rule
    cos_dtheta = np.dot(r1, r2) / (R1 * R2)
    cos_dtheta = np.clip(cos_dtheta, -1.0, 1.0)
    dtheta = np.arccos(cos_dtheta)
    if lw:
        dtheta = 2 * np.pi - dtheta
    
    if dtheta > np.pi:
        lambda_param = -np.sqrt(1.0 - c / s)
    else:
        lambda_param = np.sqrt(1.0 - c / s)

    T = np.sqrt(2 * mu / (s**3)) * tof
    
    # 4. Helper function to calculate time of flight given x
    def tof_equation(x):
        y = np.sqrt(1.0 - lambda_param**2 * (1.0 - x**2))
        if x < 1.0:
            alpha = np.acos(x)
            beta = np.acos(y)
            # Elliptic case
            return (alpha - np.sin(alpha)*np.cos(alpha) - (beta - np.sin(beta)*np.cos(beta))) / (x**3 if x != 0 else 1.0) # conceptual formula
        else:
            # Hyperbolic case
            alpha = np.acosh(x)
            beta = np.acosh(y)
            return (np.sinh(alpha)*np.cosh(alpha) - alpha - (np.sinh(beta)*np.cosh(beta) - beta)) / (x**3)

    # Refined Izzo's universal TOF evaluator to handle all x cleanly
    def get_T(x):
        y = np.sqrt(1.0 - lambda_param**2 * (1.0 - x**2))
        if x >= -1.0 and x <= 1.0:
            psi = np.arccos(x)
            phi = np.arccos(y)
            return (psi - np.sin(psi)*np.cos(psi) - (phi - np.sin(phi)*np.cos(phi))) / (x**3) if abs(x) > 1e-6 else 2/3 * (1 - lambda_param**3)
        else:
            # x > 1 (Hyperbolic)
            psi = np.arccosh(x)
            phi = np.arccosh(y)
            return (np.sinh(psi)*np.cosh(psi) - psi - (np.sinh(phi)*np.cosh(phi) - phi)) / (x**3)

    # 5. Iterative Secant Solver for x
    # Initial guesses
    x0 = 0.0
    x1 = 0.5
    
    f0 = get_T(x0) - T
    f1 = get_T(x1) - T
    
    for _ in range(max_iter):
        if abs(f1 - f0) < 1e-15:
            break
        x_next = x1 - f1 * (x1 - x0) / (f1 - f0)
        
        if abs(x_next - x1) < tol:
            x1 = x_next
            break
            
        x0, x1 = x1, x_next
        f0 = f1
        f1 = get_T(x1) - T

    # 6. Reconstruct velocities
    x = x1
    y = np.sqrt(1.0 - lambda_param**2 * (1.0 - x**2))
    
    gamma = np.sqrt(mu * s / 2.0)
    rho = (R1 - R2) / c
    sigma = np.sqrt(1.0 - rho**2)
    
    vr1 = gamma * ((lambda_param * y - x) - rho * (lambda_param * y + x)) / R1
    vt1 = gamma * sigma * (lambda_param * y + x) / R1
    
    vr2 = gamma * ((lambda_param * y - x) + rho * (lambda_param * y + x)) / R2
    vt2 = gamma * sigma * (lambda_param * y + x) / R2
    
    v1 = vr1 * ir1 + vt1 * it1
    v2 = vr2 * ir2 + vt2 * it2
    
    return v1, v2

# === Example Usage Validation ===
if __name__ == "__main__":
    # Earth gravitational parameter (km^3/s^2)
    mu_earth = 398600.4418
    
    # Low Earth Orbit to Higher Orbit Transfer vectors
    pos1 = [5000.0, 10000.0, 2100.0]     # km
    pos2 = [-14600.0, 2500.0, 6000.0]    # km
    dt = 3600.0                          # 1 hour in seconds
    
    v1, v2 = lambert_solver(pos1, pos2, dt, mu_earth)
    
    print("Departure Velocity (v1):", v1, "km/s")
    print("Arrival Velocity (v2):  ", v2, "km/s")
