import requests
import numpy as np
import pandas as pd
import json

class MonthlyData(object):
    """Monthly data object"""

    def __init__(self, lat, lon, start_year=2005, end_year=2005, database='PVGIS-SARAH', angle=0, **kwargs):
        """Constructor"""
        # Initialization :
        # Latitude :
        self._lat = lat

        # Longitude :
        self._lon = lon

        # First and last year of the data :
        self._start_year = start_year
        self._end_year = end_year

        # Default database :
        _database_option = ['PVGIS-SARAH', 'PVGIS-CMSAF', 'PVGIS-COSMO', 'PVGIS-ERA5']

        if database not in _database_option:
            print("Wrong databse")

        self._database = database

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
    def latitude(self):
        """Return the latitude of the site"""
        # Output :
        return self._lat

    @latitude.setter
    def latitude(self, lat):
        """Set the latitude of the site"""
        # Check :

        self._lat = lat

    @property
    def longitude(self):
        """Return the longitude of the site"""
        # Output :
        return self._lon

    @longitude.setter
    def longitude(self, lon):
        """Set the longitude of the site"""
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


class DailyData(object):
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
    def latitude(self):
        """Return the latitude of the site"""
        # Output :
        return self._lat

    @latitude.setter
    def latitude(self, lat):
        """Set the latitude of the site"""
        # Check :

        self._lat = lat

    @property
    def longitude(self):
        """Return the longitude of the site"""
        # Output :
        return self._lon

    @longitude.setter
    def longitude(self, lon):
        """Set the longitude of the site"""
        self._lon = lon

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
                                     self._options["global_2axis"], self._options["clearsky_2axis"])

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


if __name__ == '__main__':

    # Define a monthly data object :
    m = MonthlyData(45, 9, 2010, 2011, 'PVGIS-SARAH', 45)
    url = m.get_url()
    print(m.get_data(url))