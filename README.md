# Nowcasting

This code is a simplified version of the open source code from "[Macroeconomic Nowcasting and Forecasting with Big Data](https://www.newyorkfed.org/research/staff_reports/sr830.html)" by Brandyn Bok, Daniele Caratelli, Domenico Giannone, Argia M. Sbordone, and Andrea Tambalotti, *Staff Reports 830*, Federal Reserve Bank of New York (prepared for Volume 10 of the *Annual Review of Economics*).

**Note:** This simplified code is written so that it should be easy to follow and implement. These simplifications comes at the expense of removing data blocks (i.e. global, soft, real, labor) and AR(1) terms for each input series. **These modifications were implemented by Seth Leonard and are not associated with the Federal Reserve Bank of New York or any of its staff.**


## Download instructions

Download the code as a ZIP file by clicking the green 'Clone or download' button and selecting 'Download ZIP'.

## File and folder description

* `data/` : example US data downloaded from [FRED](https://fred.stlouisfed.org/)
* `functions/` : functions for loading data, estimating model, and updating predictions
* `example_DFM.m` : example script to estimate a dynamic factor model (DFM) for a panel of monthly and quarterly series and to create nowcasts
* `ResDFM.mat` : example DFM estimation output
* `Spec_US_example.xls` : example model specification for the US. Note that some inputs from the original code (i.e. inputs relating to blocks) are not used.


## Required software and versioning

MATLAB is required to run the code. The code was tested in MATLAB R2015b and later versions. Functionality with earlier versions of MATLAB is not guaranteed.
