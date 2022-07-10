# MSc-Implementation
Python code written during the development of the master's thesis.

This file explains what each .py or .csv file does.

FunctionsUFH.py - has all functions needed to implement the UFH model in modelUFH.py.
modelUFH.py - models the UFH. Calculates all temperatures inside the house throught the year.
PHE1.py - does all the calculations needed to model PHE1: size, heat transfer coef, temperatures. 
PHE2.py - does all the calculations needed to model PHE2: size, heat transfer coef, temperatures. 
SCOP.py - calculates SCOP for the season.
TimevsDatavsSoil.csv - contains hourly data for ambient and soil temperature in Utrecht from 2001-2010.
TimevsTemp.py - analyses TimevsDatavsSoil.csv.
soil_temp.py - calculates mean soil temperature for a certain BHE length.
interpolateAMRhot.py - calculates mass flow rate of MCHP regarding the values obtained from the MCHP model and number of reg. being used; calculates temperature of HTF leaving MCHP on the hot side based on mass flow rate, heating demand and properties of HTF.
interpolateAMRcold.py - calculates temperature of HTF leaving MCHP on the cold side based on mass flow rate, heating demand and properties of HTF.
plot_matrices.py - plots matrices obtained from MCHP model. calculates number of reg nedeed.
interpolateCOP.py - interpolates values for flow rate and frequency obtained from the MCHP model. Calculates best COP for each pair of frequency and flow rate.
qc.csv - data for cooling demand obtained from MCHP model.
qh.csv - data for heating demand obtained from MCHP model.
