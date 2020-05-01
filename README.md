# Nowcasting

This code is a simplified version of the open source code from "[Macroeconomic Nowcasting and Forecasting with Big Data](https://www.newyorkfed.org/research/staff_reports/sr830.html)" by Brandyn Bok, Daniele Caratelli, Domenico Giannone, Argia M. Sbordone, and Andrea Tambalotti, *Staff Reports 830*, Federal Reserve Bank of New York (prepared for Volume 10 of the *Annual Review of Economics*).

**Note:** This simplified code is written so that it should be easy to follow and implement for any standard mixed frequency data set. These simplifications and the flexibility of this code comes at the expense of removing data blocks (i.e. global, soft, real, labor) and AR(1) error terms. **These modifications were implemented by Seth Leonard/OttoQuant and are not associated with the Federal Reserve Bank of New York or any of its staff.**

## Using this code

This code allows replication of results obtained in the OttoQuant interface for dynamic factor models estimated by maximum likelihood. Example parameters and data are included for those without an OttoQuant subscription. 

To replicate a frequentest dynamic factor model:


1. Download the master branch of this repo in .zip format.
2. Extract the zip file and save it to a convenient location. We’ll refer to this file as “/Nowcasting-Public-master” below. 
3. Log into the OttoQuant [interface](app.ottoquant.com), select or upload data, and estimate your model.
4. Once estimation is complete, click “DOWNLOAD INPUT DATA AND PARAMETERS”
5. Extract the zip file
6. Find the folder “data” in the Nowcasting-Public repo you downloaded in step 1. Save the (extracted) contents of the zip file to this folder. /Nowcasting-Public-master/data should now contain:
  - estimation_data_no_outliers.csv
  - prediction_data_with_outliers.csv
  - low_frq_trends.csv
  - dates.csv
  - model.json 
7. In the folder “/Nowcasting-Public-master”, open run_DFM.m
8. Click “Run”

Once the model had converged to the maximum likelihood parameter estimates you will see several tables and charts. Results are stored in the structure “Res” (in your Workspace)

## Required software and versioning

MATLAB is required to run the code. The code was tested in MATLAB R2015b and later versions. Functionality with earlier versions of MATLAB is not guaranteed.
