import numpy as np
import matplotlib.pyplot as plt


class SolarDiagram(object):
    """Solar Diagram"""

    def __init__(self, lat, lon, days=[]):
        """Constructor"""
        # Check :
        if lat < -90 or lat > 90:
            raise ValueError("Wrong value for lat")

        if lon < 0 or lon > 360:
            raise ValueError("Wrong value for lon")

        if type(days) != list:
            raise TypeError("Wrong type for days")

        # Plot one days for each month of the year :
        if not days:
            days = [17, 45, 75, 105, 135, 165, 195, 225, 255, 285, 315, 345]

         # Assignment :
        self._lat = lat
        self._lon = lon
        self._days = days

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

    @property
    def days(self):
        """Return the days in the diagram"""
        # Output :
        return self._days

    @days.setter
    def days(self, days=[]):
        """Set the days in the diagram"""
        if type(days) != list:
            raise TypeError("Wrong type for days")

        # Plot one days for each month of the year :
        if not days:
            days = [17, 45, 75, 105, 135, 165, 195, 225, 255, 285, 315, 345]

        # Assignment :
        self._days = days

    def _solar_day_path(self, day, nPoints=48):
        """Return the solar path for the specified day"""
        # Check :
        if type(nPoints) != int:
            raise TypeError("Wrong type for nPoints")

        if nPoints < 12:
            raise ValueError("Wrong value for nPoints")

        # Initialization :
        _altitude = np.zeros(nPoints)
        _azimuth = np.zeros(nPoints)

        # Sunrise hour angle :
        _omega_sr = sunrise_hour_angle(self._lat, day)

        # Sunset hour angle :
        _omega_ss = sunset_hour_angle(self._lat, day)

        # Hour vector :
        _hh = np.linspace(_omega_sr, _omega_ss, num=nPoints)

        for _i in range(nPoints):
            # Solar altitude :
            _altitude[_i] = solar_altitude2(self._lat, self._lon, day, _hh[_i])

            # Solar azimuth :
            _azimuth[_i] = solar_azimuth4(self._lat, self._lon, day, _hh[_i])

        # Output :
        return _altitude, _azimuth

    def save(self, name, f_type='png'):
        """Save the figure with a certain filename"""
        # Check :
        if type(name) != str:
            raise TypeError("Wrong type for filename")

        if type(f_type) not in ('jpeg', 'png'):
            raise ValueError("Wrong type for f_type")

        # Filename :
        filename = name + '.' + f_type

        # Saving :
        plt.savefig(filename)


class CartesianDiagram(SolarDiagram):
    """Cartesian Solar Diagram class"""

    def __init__(self, lat, lon, days=[]):
        """Constructor"""
        super().__init__(lat, lon, days)

    def plot(self, nPoints=48):
        """Plot the solar path for the specified days"""
        # Initialization :
        _y_max = 0
        _x_max = 0

        # Solar paths :
        for _i in self._days:
            # Altitude and azimuth :
            _altitude, _azimuth = self._solar_day_path(_i, nPoints)

            # Max azimuth :
            _x_max = np.max([np.max(_azimuth), _x_max])

            # Max altitude :
            _y_max = np.max([np.max(_altitude), _y_max])

            # Plot :
            plt.plot(_azimuth, _altitude, color='k')

        # Labels :
        plt.xlabel('Azimuth - [°]')
        plt.ylabel('Altitude - [°]')

        # Limits :
        plt.xlim([-_x_max - 0.1*_x_max, _x_max + 0.1*_x_max])
        plt.ylim([0, _y_max + 0.1*_y_max])

        # Title :
        plt.title('Cartesian Solar Diagram - Lat {0}° Lon {1}°'.format(self._lat, self._lon))

        # Show :
        plt.grid(True, lw=2, ls='--', c='0.75')
        plt.show()


class PolarDiagram(SolarDiagram):
    """Polar Diagram class"""

    def __init__(self, lat, lon):
        """Constructor"""
        super().__init__(lat, lon)

    def plot(self, days, nPoints):
        """Plot the solar path for the specified days"""
        # Check :
        if type(days) != list:
            raise TypeError("Wrong type for days")

        # Initialize the figure :
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

        # Solar paths :
        for _i in days:
            _altitude, _azimuth = self._solar_day_path(_i, nPoints)
            print(_altitude)
            # Plot :
            ax.plot(_azimuth, _altitude)

        # Show :
        plt.show()


class TrackingSystem(object):
    """Tracking system object"""

    def __init__(self, lat, lon, n):
        """Constructor"""
        # Check :
        if lat < -90 or lat > 90:
            raise ValueError("Wrong value for lat")

        if lon < 0 or lon > 360:
            raise ValueError("Wrong value for lon")

        # Assignment :
        self._lat = lat
        self._lon = lon
        self._n = n

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

    @property
    def nDay(self):
        """Return the day of the year"""
        # Output :
        return self._n

    @nDay.setter
    def nDays(self, n):
        """Set the day of the year"""
        # Check :
        if n < 0 or n > 365:
            raise ValueError("Wrong value for n")

        # Assignment :
        self._n = n


class NS_horiz(TrackingSystem):
    """Horizontal N-S axis with horizontal tracking"""
    def __init__(self, lat, lon, n):
        """Constructor"""
        super().__init__(lat, lon, n)

    def incident_angle(self):
        """Return the incident angle"""
        # Declination :
        delta = declination(self._n)

        # Solar altitude :
        alpha = solar_altitude2(self._lat, self._lon, self._n)

        # Computation :
        theta = np.arccos(np.sqrt(np.sin(np.deg2rad(alpha))**2 +
                                  np.cos(np.deg2rad(delta))**2*np.sin(np.deg2rad(omega))**2))

        # Output :
        return theta

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
    B = (n - 1) * 360 / 365

    delta = 180/np.pi*(0.006918 - 0.399912*np.cos(np.deg2rad(B)) + 0.070257*np.sin(np.deg2rad(B)) -
                       0.006758*np.cos(2*np.deg2rad(B)) + 0.000907*np.sin(2*np.deg2rad(B)) -
                       0.002697*np.cos(3*np.deg2rad(B)) + 0.00148*np.sin(3*np.deg2rad(B)))

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
    AST_m = LST + ET/60 - 4/60*(lon_s - lon) - DS*1

    # Output :
    return AST_m

def AST2(omega):
    """Return the actual solar time from the hour angle"""
    # Check :
    if omega < -180 or omega > 180:
        raise ValueError("Wrong value for omega")

    AST_m = omega/15 + 12

    # Output :
    return AST_m

def LST(AST_m, lon, n, DS=0):
    """Return the local sun time"""
    # Check :
    if type(AST_m) not in (tuple, int, float, np.float64, np.int):
        raise TypeError("Wrong type for AST")

    if type(AST_m) == tuple and len(AST_m) != 2:
        raise ValueError("AST must have 2 elements")

    if lon < 0:
        lon = 360 + lon

    if lon > 360:
        raise ValueError("Wrong value for lon")

    if n < 0 or n > 365:
        raise ValueError("Wrong value for n")

    # Equation of time :
    ET = equation_time(n)  # [min]

    # Convert the Local solar time :
    if type(AST) == tuple:
        AST_m = AST_m[0] + AST_m[1] / 60  # [h]

    # Standard longitude :
    lon_s = 15 * np.ceil(lon / 15)

    # Computation :
    LST = AST_m - ET / 60 + 4 / 60 * (lon_s - lon) + DS * 1

    # Output :
    return AST_m

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

def solar_altitude2(lat, lon, n, omega):
    """Return the solar altitude"""
    # Check :
    if lat < -90 or lat > 90:
        raise ValueError("Wrong value for lat")

    if n <= 0 or n > 365:
        raise ValueError("Wrong value for n")

    # Declination of the sun :
    delta = declination(n)

    # Computation :
    alpha = np.rad2deg(np.arcsin(np.sin(np.deg2rad(lat)) * np.sin(np.deg2rad(delta))
                                 + np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(delta)) * np.cos(np.deg2rad(omega))))

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

def solar_azimuth4(lat, lon, n, omega):
    """Return the solar azimuth"""
    # Declination of the sun :
    delta = declination(n)

    # Solar altitude :
    alpha = solar_altitude2(lat, lon, n, omega)

    # Computation :
    gamma = np.rad2deg(
        np.sign(omega) * np.arccos((np.cos(np.deg2rad(delta)) * np.cos(np.deg2rad(omega)) * np.sin(np.deg2rad(lat)) -
                                    np.sin(np.deg2rad(delta)) * np.cos(np.deg2rad(lat))) / np.cos(np.deg2rad(alpha))))

    # Output :
    return gamma

def sunrise(lat, n, DS=0):
    """Return the sunrise time"""
    # Check :
    if lat < -90 or lat > 90:
        raise ValueError("Wrong value for lat")

    if n < 0 or n > 365:
        raise ValueError("Wrong value for n")

    # Declination :
    delta = spencer(n)

    # Computation :
    h_sr = 12 - 1/15*np.rad2deg(np.arccos(-np.tan(np.deg2rad(lat))*np.tan(np.deg2rad(delta))))

    # Daylight saving :
    h_sr = h_sr + DS

    # Output :
    return h_sr

def sunrise_hour_angle(lat, n):
    """Return the hour angle of the sunrise"""
    # Check :
    if lat < -90 or lat > 90:
        raise ValueError("Wrong value for lat")

    if n < 0 or n > 365:
        raise ValueError("Wrong value for n")

    # Declination :
    delta = spencer(n)

    # Computation :
    omega_sr = -np.rad2deg(np.arccos(-np.tan(np.deg2rad(lat)) * np.tan(np.deg2rad(delta))))

    # Output :
    return omega_sr

def sunset(lat, n, DS=0):
    """Return the sunset time"""
    # Check :
    if lat < -90 or lat > 90:
        raise ValueError("Wrong value for lat")

    if n < 0 or n > 365:
        raise ValueError("Wrong value for n")

    # Declination :
    delta = spencer(n)

    # Computation :
    h_ss = 12 + 1 / 15 * np.rad2deg(np.arccos(-np.tan(np.deg2rad(lat)) * np.tan(np.deg2rad(delta))))

    # Daylight saving :
    h_ss = h_ss + DS

    # Output :
    return h_ss

def sunset_hour_angle(lat, n):
    """Return the hour angle of the sunset"""
    # Check :
    if lat < -90 or lat > 90:
        raise ValueError("Wrong value for lat")

    if n < 0 or n > 365:
        raise ValueError("Wrong value for n")

    # Declination :
    delta = spencer(n)

    # Computation :
    omega_ss = np.rad2deg(np.arccos(-np.tan(np.deg2rad(lat)) * np.tan(np.deg2rad(delta))))

    # Output :
    return omega_ss


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
