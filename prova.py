import requests
import pandas as pd
import json

# Latitude :
lat = 45

# Longitude :
lon = 9

# Start year :
start_year = 2010

# End year :
end_year = 2011

# Database :
database = 'PVGIS-SARAH'

# Inclination :
angle = 45

# Output format :
out_format = 'json'

# Option horizontal irradiance :
option_horirrad=1

# Option direct irradiance :
option_mr_dni = 1

# Option optimal angle radiation :
option_optrad = 1

# Option global radiation at a specified angle :
option_selectrad = 1

# Option temperature :
option_temp = 1

# Base url :
base_url = 'https://re.jrc.ec.europa.eu/api/MRcalc?'

url_m = "lat={0}&lon={1}&startyear={2}&endyear={3}&raddatabase={4}&angle={5}&" \
        "browser=1&outputformat=json&select_database_month={4}&mstartyear={2}&" \
        "mendyear={3}&horirrad={6}&mr_dni={7}&optrad={8}&selectrad={9}&mangle={5}&" \
        "avtemp={10}".format(lat, lon, start_year, end_year, database, angle, option_horirrad,
                           option_mr_dni, option_optrad, option_selectrad, option_temp)

url = base_url + url_m

# Request the base site :
r = requests.get(url)

print(r.text)

result = json.loads(r.text)

print(result["outputs"]["monthly"])

class MonthlyData(object):
    """Monthly data object"""

    def __init__(self, lat, lon, start_year, end_year, database, angle, option):
        """Constructor"""
        # Initialization :
        self._lat = lat
        self._lon = lon
        self._start_year = start_year
        self._end_year = end_year
        self._database = database
        self._angle = angle
        self._option = option

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
        return self._option

    @option.setter
    def option(self, option):
        """Set the option"""
        self._option = option

    def get_url(self):
        # Base url :
        base_url = 'https://re.jrc.ec.europa.eu/api/MRcalc?'

        url_m = "lat={0}&lon={1}&startyear={2}&endyear={3}&raddatabase={4}&angle={5}&" \
                "browser=1&outputformat=json&select_database_month={4}&mstartyear={2}&" \
                "mendyear={3}&horirrad={6}&mr_dni={7}&optrad={8}&selectrad={9}&mangle={5}&" \
                "avtemp={10}".format(self._lat, self._lon, self._start_year, self._end_year,
                                     self._database, self._angle, option_horirrad,
                                     option_mr_dni, option_optrad, option_selectrad, option_temp)
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

        # Insert the data into the pandas Dataframe :
        data = pd.DataFrame()

        # Initialization :
        _years_col = []
        _month_col = []
        _H_hor_col = []
        _H_opt_col = []
        _H_i_col = []

        for i in range(0, len(result["outputs"]["monthly"])):
            _years_col.append(result["outputs"]["monthly"][i]["year"])
            _month_col.append(result["outputs"]["monthly"][i]["month"])
            _H_hor_col.append(result["outputs"]["monthly"][i]["H(h)_m"])
            _H_opt_col.append(result["outputs"]["monthly"][i]["H(i_opt)_m"])
            _H_i_col.append(result["outputs"]["monthly"][i]["H(i)_m"])

        # Output :
        return _years_col

    def __str__(self):
        """Print method"""
        pass

class DailyData(object):
    """Daily data class"""

    def __init__(self):
        """Constructor"""
        pass

if __name__ == '__main__':

    # Define a monthlydata object :
    m = MonthlyData(45, 9, 2010, 2011, 'PVGIS-SARAH', 45, 10)
    url = m.get_url()
    print(m.get_data(url))