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

print(result["outputs"]["monthly"][0]["year"])

class MonthlyData(object):
    """Monthly data object"""

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

    def get_url(self, lat, lon, start_year, end_year, database):
        # Base url :
        base_url = 'https://re.jrc.ec.europa.eu/api/MRcalc?'

        url_m = "lat={0}&lon={1}&startyear={2}&endyear={3}&raddatabase={4}&angle={5}&" \
                "browser=1&outputformat=json&select_database_month={4}&mstartyear={2}&" \
                "mendyear={3}&horirrad={6}&mr_dni={7}&optrad={8}&selectrad={9}&mangle={5}&" \
                "avtemp={10}".format(lat, lon, start_year, end_year, database, angle, option_horirrad,
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

class DailyData(object):
    """Daily data class"""

    def __init__(self):
        """Constructor"""
        pass
