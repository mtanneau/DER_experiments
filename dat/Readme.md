# Data description

These files contain hourly price, consumption, production and temperature data for the Canadian province of Ontario, in the year 2016.

This data (as well as other historical datasets) can be obtained from the [IESO website](www.ieso.ca).

## Consumption data

`load2016.csv` contains hourly consumption data.
* `OntDemand`: Ontario's total electricity consumption, in MW.h. 

## Pricing data

`price2016.csv` contains hourly pricing data.
* `HOEP`: Hourly Ontario Energy Price, in CAD per kilowatt-hour.
* `TOU`: Hourly Time-Of-Use rates for Ontario residential customers. See the [IESO website](www.ieso.ca) for more details on TOU rates.

## Production data
`prod2016.csv` contains hourly electricity production data. All outputs are in MW.h.
* `NUCLEAR`: Nuclear power plants
* `GAS`: Gas-fired power plants
* `HYDRO`: Hydroelectric power plants
* `WIND`: Wind power plants
* `SOLAR`: Solar photovoltaic power plants
* `BIOFUEL`: Biofuel power plants


## Temperature data
`temperature2016.csv` contains hourly temperature readings from a weather station in the Toronto area. All temperatures are in degrees Celsius.
Please note that some data may be missing!