# install.packages("remotes")
remotes::install_github("SRJPE/SRJPEdata")#,force=T
remotes::install_github("SRJPE/SRJPEmodel")


data(package="SRJPEdata")#shows all tables

#Data Objects
d1=SRJPEdata::weekly_efficiency
d2=SRJPEdata::weekly_juvenile_abundance_efficiency_data