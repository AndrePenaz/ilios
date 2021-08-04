import SunRadiation
import PVGIS

def clearness_index(lat, H, n):
    """Compute the clearness index"""
    # Check :
    if H <= 0:
        raise ValueError("Wrong value for H")

    # Extraterrestrial irradiation :
    H_0 = SunRadiation.extra_irradiation(lat, n)

    # Clearness index :
    k_t = H/H_0

    # Output :
    return k_t

def clearness_index2(lat, lon, n, h):
    """Compute the clearness index from the PVGIS database data"""


def liu_jordan(k_t):
    """Liu-Jordan diffuse irradiation model"""
    # Check :
    if k_t < 0 or k_t > 1:
        raise ValueError("Wrong value for k_t")

    # Correlation:
    k = 1.39 - 4.027*k_t + 5.531*k_t**2 - 3.108*k_t**3

    # Output :
    return k


def page(k_t):
    """Page's diffuse irradiation model"""
    # Check :
    if k_t < 0 or k_t > 1:
        raise ValueError("Wrong value for k_t")

    # Correlation :
    k = 1 - 1.13*k_t

    # Output :
    return k

def erbs(k_t):
    """Erbs' diffuse irradiation model"""
    # Check :
    if k_t < 0 or k_t > 1:
        raise ValueError("Wrong value for k_t")

    # Correlation :
    k = 1.317 - 3.023*k_t + 3.372*k_t**2 - 1.769*k_t**3

    # Output :
    return k

def barbaro(k_t, nCorr=0):
    """Barbaro's diffuse irradiation model"""
    # Check :
    if k_t < 0 or k_t > 1:
        raise ValueError("Wrong value for k_t")

    if nCorr == 1:
        # North Italy :
        k = 1.0492 - 1.3246*k_t
    elif nCorr == 2:
        # Central Italy :
        k = 1.0896 - 1.4797*k_t + 0.1471*k_t**2
    elif nCorr == 3:
        # South Italy :
        k = 13.9375 - 76.276*k_t + 144.3846*k_t**2 - 92.148*k_t**3
    else:
        raise ValueError("Wrong value for nCorr")

