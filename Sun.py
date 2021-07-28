import numpy as np
import matplotlib.pyplot as plt

class SolarDiagram(object):
    """Solar Diagram class"""

    def __init__(self, lat, lon):
        """Constructor"""
        # Assignment :
        self._lat = lat

        self._lon = lon

    @property
    def latitude(self):
        """Return the latitude of the site"""
        # Output :
        return self._lat

    @latitude.setter
    def latitude(self, lat):
        """Set the latitude of the site"""
        # Check :
        if lat < -90 or lat > 90:
            raise ValueError("Wrong value for lat")

        # Assignment :
        self._lat = lat

    @property
    def longitude(self):
        """Return the longitude of the site"""
        # Output :
        return self._lon

    @longitude.setter
    def longitude(self, lon):
        """Return the longitude of the site"""
        # Check :
        if lon < 0 or lon > 360:
            raise ValueError("Wrong value for lon")

        # Assignment :
        self._lon = lon

    def _solar_day_path(self, day, nPoints=24):
        """Return the solar path for the specified day"""
        # Check :
        if type(nPoints) != int:
            raise TypeError("Wrong type for nPoints")

        if nPoints < 12:
            raise ValueError("Wrong value for nPoints")

        # Initialization :
        _altitude = np.zeros(nPoints)
        _azimuth = np.zeros(nPoints)

        # Sunrise :
        _h_sr = sunrise(self._lat, day)
        _omega_sr = hour_angle(_h_sr, self._lon, day)

        # Sunset :
        _h_ss = sunset(self._lat, day)
        _omega_ss = hour_angle(_h_ss, self._lon, day)

        # Hour vector :
        _hh = np.linspace(_h_sr, _h_ss, num=nPoints)

        for _i in range(nPoints):
            # Solar altitude :
            _altitude[_i] = solar_altitude(self._lat, self._lon, day, _hh[_i])

            # Solar azimuth :
            _azimuth[_i] = solar_azimuth3(self._lat, self._lon, day, _hh[_i])

        # Output :
        return (_altitude, _azimuth)

    def plot(self, days, nPoints):
        """Plot the solar path for the specified days"""
        # Check :
        if type(days) != list:
            raise TypeError("Wrong type for days")

        # Initialize the figure :


        # Solar paths :
        for _i in days:
            _altitude, _azimuth = self._solar_day_path(_i, nPoints)
            print(_altitude)
            # Plot _
            plt.plot(_azimuth, _altitude)

        # Show :
        plt.show()


# Solar constant :
G_sc = 1352                          # [W/m^2]

def extra_radiation(n):
    """Return the extraterrestrial radiation for the specified day of the year"""
    # Check :
    if n < 0 or n > 365:
        raise ValueError("Wrong value for n")

    # Computation :
    B = (n-1)*360/365

    G = G_sc*(1.000110 + 0.034221*np.cos(np.deg2rad(B)) + 0.001280*np.sin(np.deg2rad(B)) +
              0.00719*np.cos(2*np.deg2rad(B)) + 0.000077*np.sin(2*np.deg2rad(B)))

    # Output :
    return G

def date2n(day, month):
    """Return the number of the day in the year from the data"""
    # Check :

    months_days = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

    # Computation :
    n = sum(months_days[0:month-1]) + day

    # Output :
    return n

def equation_time(n):
    """Compute the equation of time"""
    # Check :
    if n < 0 or n > 365:
        raise ValueError("Wrong value for n")

    # Computation :
    B = (n - 1) * 360 / 365

    E = 229.2*(0.000075 + 0.001868*np.cos(np.deg2rad(B)) - 0.032077*np.sin(np.deg2rad(B)) -
        0.014615*np.cos(2*np.deg2rad(B)) - 0.04089*np.sin(2*np.deg2rad(B)))

    # Output :
    return E

def declination(n):
    """Compute the declination of the sun"""
    # Check :
    if n < 0 or n > 365:
        raise ValueError("Wrong value for n")

    # Computation :
    delta = 23.45*np.sin(np.deg2rad(360*(284+n)/365))

    # Output :
    return delta

def spencer(n):
    """Compute the declination of the sun based on Spencer formula"""
    # Check :
    if n < 0 or n > 365:
        raise ValueError("Wrong value for n")

    # Computation :
    gamma = 2*np.pi*(n-1)/365

    delta = (0.006918 - 0.399912*np.cos(gamma) + 0.070257*np.sin(gamma) -
             0.006758*np.cos(2*gamma) + 0.000907*np.sin(2*gamma)  - 0.002697*np.cos(3*gamma)
             + 0.00148*np.sin(3*gamma))

    # Output :
    return delta

def AST(LST, lon, n, DS=0):
    """Return the apparent solar time for the specified location"""
    # Check :
    if type(LST) not in (tuple, int, float, np.float64, np.int):
        raise TypeError("Wrong type for LST")

    if type(LST) == tuple and len(LST) != 2:
        raise ValueError("LST must have 2 elements")

    if lon < 0:
        lon = 360 + lon

    if lon > 360:
        raise ValueError("Wrong value for lon")

    if n < 0 or n > 365:
        raise ValueError("Wrong value for n")

    # Equation of time :
    ET = equation_time(n)                    # [min]

    # Convert the Local solar time :
    if type(LST) == tuple:
        LST = LST[0] + LST[1]/60                 # [h]

    # Standard longitude :
    lon_s = 15*np.ceil(lon/15)

    # Computation :
    AST = LST + ET/60 - 4/60*(lon_s - lon) - DS*1

    # Output :
    return AST

def hour_angle(LST, lon, n, DS=0):
    """Return the hour angle at the specified location"""
    # Apparent solar time :
    AST_n = AST(LST, lon, n, DS)

    # Computation :
    omega = 15*(AST_n - 12)

    # Output :
    return omega

def solar_altitude(lat, lon, n, LST, DS=0):
    """Return the altitude of the sun"""
    # Check :
    if lat < -90 or lat > 90:
        raise ValueError("Wrong value for lat")

    if n <= 0 or n > 365:
        raise ValueError("Wrong value for n")

    # Declination of the sun :
    delta = declination(n)

    # Hour angle :
    omega = hour_angle(LST, lon, n, DS=0)

    # Computation :
    alpha = np.rad2deg(np.arcsin(np.sin(np.deg2rad(lat))*np.sin(np.deg2rad(delta))
                      + np.cos(np.deg2rad(lat))*np.cos(np.deg2rad(delta))*np.cos(np.deg2rad(omega))))

    # Condition of the Sun under the horizon :
    if alpha < 0:
        alpha = 0

    # Output :
    return alpha

def solar_zenith(lat, lon, n, LST, DS=0):
    """Return the zenith angle of the sun"""
    # Computation :
    alpha = solar_altitude(lat,lon, n, LST, DS=0)
    
    phi = 90 - alpha

    # Output :
    return phi

def solar_azimuth(lat, lon, n, LST):
    """Return the azimuth of the sun"""
    # Declination of the sun :
    delta = declination(n)

    print("Declination :")
    print(delta)

    # Hour angle :
    omega = hour_angle(LST, lon, n, DS=0)

    print("Omega :")
    print(omega)

    # Solar altitude :
    alpha = solar_altitude(lat, lon, n, LST)

    print("Altitude")
    print(alpha)

    # Computation :
    gamma = np.rad2deg(np.arcsin(np.cos(np.deg2rad(delta))*np.sin(np.deg2rad(omega))/np.cos(np.deg2rad(alpha))))

    if np.cos(np.deg2rad(omega)) > np.tan(np.deg2rad(delta))/np.tan(np.deg2rad(lat)):
        if omega < 0:
            gamma = np.abs(gamma) - 180
        else:
            gamma = 180 - gamma

    # Output :
    return gamma

def solar_azimuth2(lat, lon, n, LST):
    delta = np.deg2rad(declination(n))
    omega = np.deg2rad(hour_angle(LST, lon, n))
    zenith = np.deg2rad(solar_zenith(lat,lon, n, LST))

    gamma = np.sign(omega)*np.abs(np.arccos((np.cos(zenith)*np.sin(np.deg2rad(lat)) - np.sin(delta))/np.sin(zenith)*np.cos(np.deg2rad(lat))))

    print(np.rad2deg(gamma))

def solar_azimuth3(lat, lon, n, LST):
    """Return the azimuth of the sun"""
    # Declination of the sun :
    delta = declination(n)

    # Hour angle :
    omega = hour_angle(LST, lon, n, DS=0)

    # Solar altitude :
    alpha = solar_altitude(lat, lon, n, LST)

    # Computation :
    gamma = np.rad2deg(np.sign(omega)*np.arccos((np.cos(np.deg2rad(delta))*np.cos(np.deg2rad(omega))*np.sin(np.deg2rad(lat)) -
              np.sin(np.deg2rad(delta))*np.cos(np.deg2rad(lat)))/np.cos(np.deg2rad(alpha))))

    # Output :
    return gamma

def sunrise(lat, n):
    """Return the sunrise time"""
    # Check :
    if lat < -90 or lat > 90:
        raise ValueError("Wrong value for lat")

    if n < 0 or n > 365:
        raise ValueError("Wrong value for n")

    # Declination :
    delta = declination(n)

    # Computation :
    h_sr = 12 - 1/15*np.rad2deg(np.arccos(-np.tan(np.deg2rad(lat))*np.tan(np.deg2rad(delta))))

    # Output :
    return h_sr

def sunset(lat, n):
    """Return the sunset time"""
    # Check :
    if lat < -90 or lat > 90:
        raise ValueError("Wrong value for lat")

    if n < 0 or n > 365:
        raise ValueError("Wrong value for n")

    # Declination :
    delta = declination(n)

    # Computation :
    h_ss = 12 + 1 / 15 * np.rad2deg(np.arccos(-np.tan(np.deg2rad(lat)) * np.tan(np.deg2rad(delta))))

    # Output :
    return h_ss

def day_length(lat, n):
    """Return the length of the day in a specified location"""
    # Check :
    if lat < -90 or lat > 90:
        raise ValueError("Wrong value for lat")

    if n < 0 or n > 365:
        raise ValueError("Wrong value for n")

    # Declination :
    delta = declination(n)

    # Computation :
    h_d = 2/15*np.rad2deg(np.arccos(-np.tan(np.deg2rad(lat)) * np.tan(np.deg2rad(delta))))

    # Output :
    return h_d

def incident_angle(lat, lon, n, LST, slope, azimuth_s):
    """Return the incident angle for the sun beam at a specified location
       and time"""
    # Check
    if slope < 0 or slope > 90:
        raise ValueError("Wrong value for slope")

    if azimuth_s < -180 or azimuth_s > 180:
        raise ValueError("Wrong value for azimuth")

    # Solar declination :
    delta = declination(n)

    # Solar altitude :
    alpha = solar_altitude(lat, lon, n, LST)

    # Hour angle :
    omega = hour_angle(LST, lon, n)

    # Incidence angle :
    theta = np.rad2deg(np.arccos(np.sin(np.deg2rad(lat))*np.sin(np.deg2rad(delta))*np.cos(np.deg2rad(slope)) -
                      np.cos(np.deg2rad(lat))*np.sin(np.deg2rad(delta))*np.sin(np.deg2rad(slope))*np.cos(np.deg2rad(azimuth_s)) +
                      np.cos(np.deg2rad(lat))*np.cos(np.deg2rad(delta))*np.cos(np.deg2rad(omega))*np.cos(np.deg2rad(slope)) +
                      np.sin(np.deg2rad(lat))*np.cos(np.deg2rad(delta))*np.cos(np.deg2rad(omega))*np.sin(np.deg2rad(slope))*np.cos(np.deg2rad(azimuth_s)) +
                      np.cos(np.deg2rad(delta))*np.sin(np.deg2rad(omega))*np.sin(np.deg2rad(slope))*np.sin(np.deg2rad(azimuth_s))))

    print("Latitude")
    print(lat)

    print("Declination :")
    print(delta)

    print("Slope :")
    print(slope)

    print("Azimuth :")
    print(azimuth_s)

    print("Omega :")
    print(omega)

    # Output :
    return theta