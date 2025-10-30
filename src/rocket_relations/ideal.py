import numpy as np

def characteristic_velocity(gamma, R, T0):
    """
    Calculate the ideal characteristic velocity (c*) for a rocket.

    Parameters
    ----------
    gamma: Ratio of specific heats (must be > 1, typically < 1.8).
    R: Specific gas constant [J/(kg·K)] (must be > 0).
    T0: Stagnation temperature [K] (must be > 0).

    Returns
    -------
    c_star: Ideal characteristic velocity [m/s].

    Raises
    ------
    TypeError
        If any input is not numeric.
    ValueError
        If gamma ≤ 1, R ≤ 0, or T0 ≤ 0.
    """
    # Convert to numpy arrays for elementwise ops
    gamma = np.asarray(gamma)
    R = np.asarray(R)
    T0 = np.asarray(T0)

    # Input validation
    if not (np.issubdtype(gamma.dtype, np.number) and 
            np.issubdtype(R.dtype, np.number) and 
            np.issubdtype(T0.dtype, np.number)):
        raise TypeError("All inputs must be numeric.")

    if np.any(gamma <= 1):
        raise ValueError("gamma must be > 1.")
    if np.any(R <= 0):
        raise ValueError("Specific gas constant R must be > 0.")
    if np.any(T0 <= 0):
        raise ValueError("Stagnation temperature T0 must be > 0.")

    # Compute c*
    c_star = (1 / gamma) * np.sqrt(
        (2 / (gamma + 1)) * ((gamma + 1) / 2) ** ((gamma + 1) / (gamma - 1)) * R * T0
    )

    return c_star

    def thrust_coefficient(gamma, pe_p0, pa_p0, Ae_Astar):
    """
    Calculate the ideal thrust coefficient (CF) for a rocket nozzle.

    Parameters
    ----------
    gamma: Ratio of specific heats (must be > 1, typically < 1.8).
    pe_p0: Exit-to-chamber pressure ratio (pe/p0), must be in [0, 1).
    pa_p0: Ambient-to-chamber pressure ratio (pa/p0), must be in [0, 1).
    Ae_Astar: Nozzle area ratio (Ae/A*), must be ≥ 1.

    Returns
    -------
    CF: Ideal thrust coefficient (dimensionless).

    Raises
    ------
    TypeError
        If any input is not numeric.
    ValueError
        If inputs violate physical domain constraints.
    """
    gamma = np.asarray(gamma)
    pe_p0 = np.asarray(pe_p0)
    pa_p0 = np.asarray(pa_p0)
    Ae_Astar = np.asarray(Ae_Astar)

    # Input validation
    if not all(np.issubdtype(x.dtype, np.number) for x in [gamma, pe_p0, pa_p0, Ae_Astar]):
        raise TypeError("All inputs must be numeric.")

    if np.any(gamma <= 1):
        raise ValueError("gamma must be > 1.")
    if np.any((pe_p0 < 0) | (pe_p0 >= 1)):
        raise ValueError("pe/p0 must be in the range [0, 1).")
    if np.any((pa_p0 < 0) | (pa_p0 >= 1)):
        raise ValueError("pa/p0 must be in the range [0, 1).")
    if np.any(Ae_Astar < 1):
        raise ValueError("Ae/A* must be ≥ 1.")

    # Compute CF
    term1 = np.sqrt(
        (2 * gamma**2 / (gamma - 1))
        * (2 / (gamma + 1)) ** ((gamma + 1) / (gamma - 1))
        * (1 - pe_p0 ** ((gamma - 1) / gamma))
    )
    term2 = (pe_p0 - pa_p0) * Ae_Astar

    CF = term1 + term2
    return CF

    
