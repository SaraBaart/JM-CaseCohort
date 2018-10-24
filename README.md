# JM-CaseCohort
Fitting joint models in a case-cohort design

This repository contains code that can be used to fit joint models for longitudinal and time-to-event data in a case-cohort design. 

The code assumes that the data is provided in long format, where the longitudinal measurements for the patients outside the case-cohort design are put to NA. The times on which these measurements are taken, however are still inlcuded in the data set. The data set therefore has the same number of rows as the hypothetical data set if the full cohort was measured. 

The models are estimated in the Bayesian framework using JAGS, thefore this software need to be installed.

'CaseCohort' contains the main code.
