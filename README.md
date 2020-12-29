# Cancer Signature Project

## *aapc_past folder*
**Objective: Calculating Average Annual Percent Change (AAPC) using data directly abstracted from SEER database** 
- **aapc_past_processing.R:** This file processes SEER excel file and gives an output ready to be used by Joinpoint Regression Program.
- **AAPC_calculation_incidence.jpt:** This is a template file used in the AAPC calculation by Joinpoint Regression Program for cancer incidence rate.
- **AAPC_calculation_mortality.jpt:** This is a template file used in the AAPC calculation by Joinpoint Regression Program for cancer mortality rate.

## *apc folder*
**Objective: Forecasting cancer incidence/mortality rates using Age-Period-Cohort model** 
- **apc_input_processing.R:** This file processes SEER excel file and gives an output ready to be used by APC webtool (https://analysistools.cancer.gov/apc/).
- **apc_forecast.R:** This file performs forecasting of cancer incidence/mortality rates using Age-Period-Cohort model. We specifically followed the method detailed in Rosenberg et al. (https://pubmed.ncbi.nlm.nih.gov/26063794/). It should be noted also that this file also contained the code written by Dr. Wenjiang J. Fu, which utilizes AutoRegressive Integrated Moving Average (ARIMA) projection to perform forecasting. 
- **apc_postprocessing.R:** This file processes the output from apc_forecast.R and gives the result with cancer incidence/mortality rates and corresponding confidence intervals.

## *aapc_future folder*
**Objective: Calculating Average Annual Percent Change (AAPC) using past data directly abstracted from SEER database and future data obtained from age-period-cohort models** 
- **aapc_calculation.R:** This file calculates AAPC using the output form apc_forecast.R. We utilize *segmented* package to calculate AAPC as Joinpoint Regression Program is not executable on Linux platform (unable to perform parallel computing). 
