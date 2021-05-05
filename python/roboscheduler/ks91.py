import numpy as np


def X_Z(Z):
    """a function to calculate equation 3 of KS91
       X gives the optical pathlength through the atmosphere at a given
       zenith distance Z (where Z is given in degrees)
    """
    X = np.power(1 - 0.96*np.sin(Z*np.pi/180.)*np.sin(Z*np.pi/180.), -0.5)

    return X


def B_zee(B_zen, kappa, Z):
    """a function to calculate equation 2 of KS91, which gives
       the dark time sky brightness as a function of:

    Parameters:
        ----------

        B_zen: np.float64
            the zenith sky brightness (in nanoLamberts)

        kappa: np.float64
            extinction coefficient (magnitudes per airmass)

        Z: np.float64
            zenith distance of the sight line (in degrees)
    """

    # calculate X_Z via a function call
    X = X_Z(Z)

    # calculate equation 2 of KS91
    B_0 = X*B_zen*np.power(10, -0.4*kappa*(X-1))

    return B_0


def KS91_deltaV(alpha, rho, zee, zee_m, B_zen=79, kappa=0.172):
    """Computer KS91 delta V for a field caused by the moon

    Parameters:
    ----------
    alpha : float, np.float64
        lunar phase angle (degress, 0 is full, 180 is new)
    rho : float, np.float64 (can be array)
        moon/sky seperation (degrees)
    zee : float, np.float64 (can be array)
        zenith distance of the sight line (degrees)
    zee_m : float, np.float64
        zenith distance of the moon (degrees)
    B_zen : float, np.float64
        the zenith sky brightness (in nanoLamberts)
    kappa : float, np.float64
        extinction coefficient (magnitudes per airmass)

    Returns:
    -------
    deltaV : np.float64 (or array)
        change in V mag at targ location due to moon
    """
    # calculate Istar (Eq. 20 in KS91)
    Istar = 10.**(-0.4*(3.84+0.026*np.abs(alpha)+(4e-9)*(alpha**4)))

    # calculate f of rho (eq. 21 in KS91)
    f_rho = (10.**5.36) * (1.06 + (np.cos(rho*(3.1415/180.)))**2.) \
              + np.power(10,  6.15-(rho/40.))

    # calculate B_moon
    B_moon = f_rho * Istar * (10.**(-0.4*kappa*X_Z(zee_m)))\
                   * (1 - np.power(10, -0.4*kappa*X_Z(zee)))

    # calculate deltaV
    deltaV = -2.5 * np.log10((B_moon + B_zee(B_zen, kappa, zee))/B_zee(B_zen, kappa, zee))

    return deltaV
