import requests
import numpy as np
import pandas as pd
import json

class PVGISData(object):
    """PVGIS data class"""

    def __init__(self):
        """Constructor"""
        pass

    @property
    def latitude(self):
        """Return the latitude of the site"""
        # Output :
        return self._lat

    @latitude.setter
    def latitude(self, lat):
        """Set the latitude of the site"""
        # Check :
        if lat > 90 or lat < -90:
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
        """Set the longitude of the site"""
        # Check :
        if lon < 0 or lon > 360:
            raise ValueError("Worng value for lon")

        # Assignment :
        self._lon = lon


class MonthlyData(PVGISData):
    """Monthly data object"""

    def __init__(self, lat, lon, start_year=2005, end_year=2005, database='PVGIS-SARAH', angle=0, **kwargs):
        """Constructor"""
        # Initialization :
        # Latitude :
        self._lat = lat

        # Longitude :
        self._lon = lon

        # Elevation :
        self._elevation = np.NaN

        # First and last year of the data :
        self._start_year = start_year
        self._end_year = end_year

        # Default database :
        _database_option = ['PVGIS-SARAH', 'PVGIS-CMSAF', 'PVGIS-COSMO', 'PVGIS-ERA5']

        if database not in _database_option:
            print("Wrong databse")

        # Assignment :
        self._database = database

        # Angle :
        self._angle = angle

        # Default options :
        _options = {'hor_irrad': 1, 'mr_dni': 1, 'opt_rad': 1, 'select_rad': 1, 'temp': 1}

        # Difference :
        _diff = set(kwargs.keys()) - set(_options.keys())
        if _diff:
            print("Wrong input option")

        # Update the options' dictionary :
        _options.update(kwargs)
        self._options = _options

    @property
    def elevation(self):
        """Return the elevation of the site"""
        self._elevation

    @property
    def start_year(self):
        """Return the start year for the data"""
        # Output :
        return self._start_year

    @start_year.setter
    def start_year(self, start_year):
        """Set the start year for the data"""
        self._start_year = start_year

    @property
    def end_year(self):
        """Return the end year for the data"""
        # Output :
        return self._end_year

    @end_year.setter
    def end_year(self, end_year):
        """Set the end year for the data"""
        self._end_year = end_year

    @property
    def database(self):
        """Return the database"""
        # Output :
        return self._database

    @database.setter
    def database(self, database):
        """Set the database"""
        # Default database :
        _database_option = ['PVGIS-SARAH', 'PVGIS-CMSAF', 'PVGIS-COSMO', 'PVGIS-ERA5']

        # Check :
        if database not in _database_option:
            print("Wrong databse")

        # Assignment :
        self._database = database

    @property
    def angle(self):
        """Return the angle value"""
        # Output :
        return self._angle

    @angle.setter
    def angle(self, angle):
        """Set the angle value"""
        self._angle = angle

    @property
    def option(self):
        """Return the option"""
        # Output :
        return self._options

    @option.setter
    def option(self, **kwargs):
        """Set the option"""
        # Difference :
        diff = set(kwargs.keys()) - set(self._options.keys())
        if diff:
            print("Wrong input option")

        # Update the options' dictionary :
        self._options.update(kwargs)

    def get_url(self):
        # Base url :
        base_url = 'https://re.jrc.ec.europa.eu/api/MRcalc?'

        url_m = "lat={0}&lon={1}&startyear={2}&endyear={3}&raddatabase={4}&angle={5}&" \
                "browser=1&outputformat=json&select_database_month={4}&mstartyear={2}&" \
                "mendyear={3}&horirrad={6}&mr_dni={7}&optrad={8}&selectrad={9}&mangle={5}&" \
                "avtemp={10}".format(self._lat, self._lon, self._start_year, self._end_year,
                                     self._database, self._angle, self._options["hor_irrad"],
                                     self._options["mr_dni"], self._options["opt_rad"],
                                     self._options["select_rad"], self._options["temp"])
        # Url :
        url = base_url + url_m

        # Output :
        return url

    def get_data(self, url):
        """Return the data"""
        # Request the base site :
        r = requests.get(url)

        # Json file :
        result = json.loads(r.text)

        print(r.text)

        # Elevation :
        self._elevation = result["inputs"]["location"]["elevation"]

        # Number of rows :
        _nRows = len(result["outputs"]["monthly"])

        # Initialization :
        _data = np.zeros((_nRows, 7))

        for i in range(0, len(result["outputs"]["monthly"])):
            _data[i][0] = result["outputs"]["monthly"][i]["year"]
            _data[i][1] = result["outputs"]["monthly"][i]["month"]
            _data[i][2] = result["outputs"]["monthly"][i]["H(h)_m"]
            _data[i][3] = result["outputs"]["monthly"][i]["H(i_opt)_m"]
            _data[i][4] = result["outputs"]["monthly"][i]["H(i)_m"]
            _data[i][5] = result["outputs"]["monthly"][i]["Hb(n)_m"]
            _data[i][6] = result["outputs"]["monthly"][i]["T2m"]

        # Insert the data into the pandas Dataframe :
        data = pd.DataFrame(_data, columns=["year", "month", "H_hor", "H_opt", "H(i)", "Hb", "T_m"])

        # Output :
        return data

    def __str__(self):
        """Print method"""
        pass


class DailyData(PVGISData):
    """Daily data class"""

    def __init__(self, lat, lon, month, slope, azimuth, database, **kwargs):
        """Constructor"""
        # Initialization :
        # Latitude :
        self._lat = lat

        # Longitude :
        self._lon = lon

        # Month :
        self._month = month

        # Slope :
        self._slope = slope

        # Azimuth :
        self._azimuth = azimuth

        # Database :
        self._database = database

        # Default options :
        _options = {'localtime': 0, 'global': 1, 'clearsky': 1, 'glob_2axis': 1, 'clearsky_2axis': 1, 'showtemperatures': 1}

        # Difference :
        _diff = set(kwargs.keys()) - set(_options.keys())
        if _diff:
            print("Wrong input option")

        # Update the options' dictionary :
        _options.update(kwargs)
        self._options = _options

    @property
    def month(self):
        """Return the month"""
        # Output :
        return self._month

    @property
    def database(self):
        """Return the database"""
        # Output :
        return self._database

    @database.setter
    def database(self, database):
        """Set the database"""
        # Default database :
        _database_option = ['PVGIS-SARAH', 'PVGIS-CMSAF', 'PVGIS-COSMO', 'PVGIS-ERA5']

        if database not in _database_option:
            print("Wrong databse")

        # Assignment :
        self._database = database

    @month.setter
    def month(self, month):
        """Set the month"""
        # Assignment :
        self._month = month

    @property
    def slope(self):
        """Return the slope of the surface"""
        # Output :
        return self._slope

    @slope.setter
    def slope(self, slope):
        """Set the slope of the surface"""
        # Assignment :
        self._slope = slope

    @property
    def azimuth(self):
        """Return the azimuth of the surface"""
        # Output :
        return self._azimuth

    @azimuth.setter
    def azimuth(self, azimuth):
        """Set the azimuth of the surface"""
        # Assignment :
        self._azimuth = azimuth

    @property
    def option(self):
        """Return the option"""
        # Output :
        return self._options

    @option.setter
    def option(self, **kwargs):
        """Set the option"""
        # Difference :
        diff = set(kwargs.keys()) - set(self._options.keys())
        if diff:
            print("Wrong input option")

        # Update the options' dictionary :
        self._options.update(kwargs)

    def get_url(self):
        # Base url :
        base_url = 'https://re.jrc.ec.europa.eu/api/DRcalc?'

        url_m = "lat={0}&lon={1}&raddatabase={2}&angle={3}&aspect{4}&browser=1&outputformat=json&" \
                "select_database_daily={2}&month={5}&localtime={6}&global={7}&" \
                "clearsky={8}&dangle={3}&daspect={4}&global_2axis={8}&clearsky_2axis={9}&" \
                "showtemperatures={10}".format(self._lat, self._lon, self._database, self._slope, self._azimuth,
                                     self._month, self._options["localtime"], self._options["global"], self._options["clearsky"],
                                     self._options["glob_2axis"], self._options["clearsky_2axis"])

        # Url :
        url = base_url + url_m

        # Output :
        return url

    def get_data(self, url):
        """Return the data"""
        # Request the base site :
        r = requests.get(url)

        # Json file :
        result = json.loads(r.text)

        print(r.text)

        # Elevation :
        self._elevation = result["inputs"]["location"]["elevation"]

        # Number of rows :
        _nRows = len(result["outputs"]["daily_profile"])

        # Initialization :
        _month_col, _time_col, _G_i, _Gb, _Gd, _Gcs_i, _Gcs_n, _T2m = ([] for i in range(8))

        for i in range(0, _nRows):
            _month_col.append(result["outputs"]["daily_profile"][i]["month"])
            _time_col.append(result["outputs"]["daily_profile"][i]["time"])
            _G_i.append(result["outputs"]["daily_profile"][i]["G(i)"])
            _Gb.append(result["outputs"]["daily_profile"][i]["Gb(i)"])
            _Gd.append(result["outputs"]["daily_profile"][i]["Gd(i)"])
            _Gcs_i.append(result["outputs"]["daily_profile"][i]["Gcs(i)"])
            _Gcs_n.append(result["outputs"]["daily_profile"][i]["Gcs(n)"])
            _T2m.append(result["outputs"]["daily_profile"][i]["T2m"])

        _data = {"month": _month_col, "time": _time_col, "G_i": _G_i, "Gb": _Gb,
                 "Gd": _Gd, "Gcs_i": _Gcs_i, "Gcs_n": _Gcs_n, "Tm": _T2m}

        # Insert the data into the pandas Dataframe :
        data = pd.DataFrame(_data)

        # Output :
        return data


class HourlyData(PVGISData):
    """Hourly data class"""

    def __init__(self, lat, lon, slope=0, azimuth=0, start_year=2005, end_year=2005, database='PVGIS-SARAH',
                 power=1, loss=14, **kwargs):
        """Constructor"""
        # Latitude and longitude :
        self._lat, self._lon = lat, lon

        # Default database :
        _database_option = ['PVGIS-SARAH', 'PVGIS-CMSAF', 'PVGIS-COSMO', 'PVGIS-ERA5']

        if database not in _database_option:
            print("Wrong databse")

        # Assignment :
        self._database = database

        # Slope and azimuth of the surface :
        self._slope, self._azimuth = slope, azimuth

        # First and last year :
        self._start_year = start_year
        self._end_year = end_year

        # PV peak power :
        self._power = power

        # PV percentage loss :
        self._loss = loss

        # Default options :
        _options = {'tracking': 0}

        # Difference :
        _diff = set(kwargs.keys()) - set(_options.keys())
        if _diff:
            print("Wrong input option")

        # Update the options' dictionary :
        _options.update(kwargs)
        self._options = _options

    @property
    def slope(self):
        """Return the slope of the surface"""
        # Output :
        return self._slope

    @slope.setter
    def slope(self, slope):
        """Set the slope of the surface"""
        # Check :
        if slope < 0 or slope > 90:
            raise ValueError("Wrong value for slope")

        # Assignment :
        self._slope = slope

    @property
    def azimuth(self):
        """Return the azimuth of the surface"""
        # Output :
        return self._azimuth

    @azimuth.setter
    def azimuth(self, azimuth):
        """Set the azimuth of the surface"""
        # Check :
        if azimuth < -180 or azimuth > 180:
            raise ValueError("Wrong value for azimuth")

        # Assignment :
        self._azimuth = azimuth

    @property
    def power(self):
        """Return the power of the PV plant"""
        # Output :
        return self._power

    @power.setter
    def power(self, power):
        """Set the power of the PV plant"""
        # Check :
        if power <= 0:
            raise ValueError("Wrong value for power")

        # Assignment :
        self._power = power

    @property
    def loss(self):
        """Return the percentage of the PV loss"""
        # Output :
        return self._loss

    @loss.setter
    def loss(self, loss):
        """Set the percentage of the PV loss"""
        # Check :
        if loss < 0 or loss > 100:
            raise ValueError("Wrong value for loss")

        # Assignment :
        self._loss = loss

    def get_url(self):
        """"Return the url"""
        # Base url :
        base_url = 'https://re.jrc.ec.europa.eu/api/seriescalc?'

        url_m = "lat={0}&lon={1}&raddatabase={2}&angle={3}&aspect{4}&startyear={5}&endyear={6}&" \
                "browser=1&outputformat=json&usehorizon=1&mountingplace=free&optimalinclination=1&" \
                "optimalangles=0&select_database_hourly={2}&hstartyear={5}&hendyear={6}&" \
                "trackingtype={9}&hourlyangle={3}&hourlyoptimalinclination=1&hourlyaspect={4}&" \
                "pvcalculation=1&pvtechchoice=crystSi&peakpower={7}&loss={8}&" \
                "components=1&".format(self._lat, self._lon, self._database, self._slope, self._azimuth,
                                       self._start_year, self._end_year, self._power, self._loss,
                                       self._options["tracking"])

        # Url :
        url = base_url + url_m

        # Output :
        return url

    def get_data(self, url):
        """Return the data"""
        # Request the base site :
        r = requests.get(url)

        # Json file :
        result = json.loads(r.text)

        print(r.text)

        # Elevation :
        self._elevation = result["inputs"]["location"]["elevation"]

        # Number of rows :
        _nRows = len(result["outputs"]["hourly"])

        # Initialization :
        _time, _P, _Gb, _Gd, _Gr, _H_sun, _T2m, _WS10m, _Int = ([] for i in range(9))

        for i in range(0, _nRows):
            _time.append(result["outputs"]["hourly"][i]["time"])
            _P.append(result["outputs"]["hourly"][i]["P"])
            _Gb.append(result["outputs"]["hourly"][i]["Gb(i)"])
            _Gd.append(result["outputs"]["hourly"][i]["Gd(i)"])
            _Gr.append(result["outputs"]["hourly"][i]["Gr(i)"])
            _H_sun.append(result["outputs"]["hourly"][i]["H_sun"])
            _T2m.append(result["outputs"]["hourly"][i]["T2m"])
            _WS10m.append(result["outputs"]["hourly"][i]["WS10m"])
            _Int.append(result["outputs"]["hourly"][i]["Int"])

        _data = {"time": _time, "P": _P, "Gb": _Gb, "Gd": _Gd, "Gr": _Gr,
                 "H_sun": _H_sun, "Tm": _T2m, "WS10m": _WS10m, "Int": _Int}

        # Insert the data into the pandas Dataframe :
        data = pd.DataFrame(_data)

        # Output :
        return data


class TMY(object):
    """Typical Meteorogical Year class"""

    def __init__(self, lat, lon, start_year=2005, end_year=2014):
        """Constructor"""
        # Initialization :
        # Latitude :
        self._lat = lat

        # Longitude :
        self._lon = lon

        # Elevation :
        self._elevation = np.NaN

        # First and last year of the data :
        self._start_year = start_year
        self._end_year = end_year

        # Horizon :
        self._horizon = 1

    @property
    def latitude(self):
        """Return the latitude of the site"""
        # Output :
        return self._lat

    @latitude.setter
    def latitude(self, lat):
        """Set the latitude of the site"""
        # Check :
        if lat > 90 or lat < -90:
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
        """Set the longitude of the site"""
        # Check :
        if lon < 0 or lon > 360:
            raise ValueError("Worng value for lon")

        # Assignment :
        self._lon = lon

    @property
    def start_year(self):
        """Return the start year for the data"""
        # Output :
        return self._start_year

    @start_year.setter
    def start_year(self, start_year):
        """Set the start year for the data"""
        self._start_year = start_year

    @property
    def end_year(self):
        """Return the end year for the data"""
        # Output :
        return self._end_year

    @end_year.setter
    def end_year(self, end_year):
        """Set the end year for the data"""
        self._end_year = end_year

    def get_url(self):
        """Return the url"""
        # Base url :
        base_url = 'https://re.jrc.ec.europa.eu/api/tmy?'

        url_m = 'lat={0}&lon={1}&usehorizon={2}&browser=1&outputformat=json&startyear={3}&' \
                'endyear={4}&period=2'.format(self._lat, self._lon, self._horizon, self._start_year, self._end_year)

        url = base_url + url_m

        # Output :
        return url

    def get_data(self, url):
        """Return the data"""
        # Request the base site :
        r = requests.get(url)

        # Json file :
        result = json.loads(r.text)

        print(r.text)

        # Elevation :
        self._elevation = result["inputs"]["location"]["elevation"]

        # Number of rows :
        _nRows = len(result["outputs"]["tmy_hourly"])

        # Initialization :
        _time, _T2m, _RH, _G, _Gb, _Gd, _IR, _ws10m, _wd10m, _sap  = ([] for i in range(10))

        for i in range(0, _nRows):
            _time.append(result["outputs"]["tmy_hourly"][i]["time(UTC)"])
            _T2m.append(result["outputs"]["tmy_hourly"][i]["T2m"])
            _RH.append(result["outputs"]["tmy_hourly"][i]["RH"])
            _G.append(result["outputs"]["tmy_hourly"][i]["G(h)"])
            _Gb.append(result["outputs"]["tmy_hourly"][i]["Gb(n)"])
            _Gd.append(result["outputs"]["tmy_hourly"][i]["Gd(h)"])
            _IR.append(result["outputs"]["tmy_hourly"][i]["IR(h)"])
            _ws10m.append(result["outputs"]["tmy_hourly"][i]["WS10m"])
            _wd10m.append(result["outputs"]["tmy_hourly"][i]["WD10m"])
            _sap.append(result["outputs"]["tmy_hourly"][i]["SP"])

        _data = {"time": _time, "Tm": _T2m, "RH": _RH, "G": _G,
                 "Gb": _Gb, "Gd": _Gd, "IR": _IR, "ws10m": _ws10m, "wd10m": _wd10m,
                 "sap": _sap}

        # Insert the data into the pandas Dataframe :
        data = pd.DataFrame(_data)

        # Output :
        return data


if __name__ == '__main__':

    # Define a monthly data object :
    m = MonthlyData(45, 9, 2010, 2011, 'PVGIS-SARAH', 45)
    url = m.get_url()
    print(m.get_data(url))

    # Define a daily data object :
    d = DailyData(45, 9, 6, 30, 0, 'PVGIS-SARAH')
    url = d.get_url()
    print(d.get_data(url))

    # Define a TMY data object :
    tmy = TMY(45, 9)
    url = tmy.get_url()
    print(tmy.get_data(url))

    # Define a hourly data object :
    h = HourlyData(45, 9, 45, 0)
    url = h.get_url()
    print(h.get_data(url))

