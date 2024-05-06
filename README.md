# SRJPEmodel

This data package contains the Spring Run Juvenile Production Estimate (SRJPE) Forecast model. This model is made up of multiple submodels that are combined to produce an annual forecast of juvenile spring run chinook salmon at the delta entry. 

This model is dependent on data produced in the `SRJPEdata` R package. 

The SR JPE model is still in development. This package will be updated to contain the most up to date code and documentation on the model. 

## Installation

```
# install.packages("remotes")
remotes::install_github("SRJPE/SRJPEmodel")
```

## Usage
This package contains code and documenation for the SR JPE model. 

```
# datasets within the package
data(package = 'SRJPEmodel')

# explore package documentation 
?SRJPEmodel
```

## SR JPE Submodel Components 
This repository contains code to run a few versions of the SR JPE including: 

* Within Season Forecast 
* SR Annual Forecast 

*Both these versions can be run using tributary RST data or Mainstem RST data*

These models are both dependent on the following submodels, click on the link to view an article giving a brief overview of each submodel: 

* BTSPAS-X Juvenile Abundance Modeling 
* Passage to Spawner (P2S) - Adult Spawner Modeling 
* Survival Model 
* Probabilistic Length at Date (PLAD) - Run Identification Model 
* Stock-Recruit Model 
* Travel Time Model 

## Explore model fits and results 

Explore model datasets, fits, and results on the SRJPEdashboard shiny app. (Add link)



