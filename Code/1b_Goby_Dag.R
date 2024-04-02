##DAGs for TW Goby

library(dagitty)

## Dependent Var
# TW Goby density (area or volume)
# TW Goby size distribution?

## Independent covariates
# Year
# Season (only a few early dates so probably skip)
# Water temp
# Salinity
# DO (perc_sat)
# pH
# Stickleback density
# Sculpin Density
# microsporidia frequency
# substrate
# phytoplankton cover?  May need to get into dataset from raw data
# SAV species?
# SAV Cover?
# EAV Species?
# EAV Cover?
# Zone ???  WQ and substrate act as controls, esp. low DO
# Breach <- Need to code  #DF suggests using date since last Breach...can provide dates
# Wind/Weather <- Need to get

## Random Effects
# Zone/Station_ID

# build the DAG

set.seed(123)

DAG_GOBY <- dagitty("dag{ 
  Year -> GOBY ;
  DO_Temp -> GOBY ;
  Stickleback -> GOBY ;
  Sculpin -> GOBY ;
  Microsporidia -> GOBY ;
  Substrate -> GOBY ;
  SAV -> GOBY ;
  Rain -> Breach -> DO_Temp -> GOBY ;
  Rain -> DO_Temp -> GOBY ;
  Wind -> DO_Temp -> GOBY ;
  Zone -> DO_Temp -> GOBY ;
  DO_Temp -> SAV -> GOBY ;
  DO_Temp -> Sculpin -> GOBY ;
  DO_Temp -> Stickleback -> GOBY ;
  DO_Temp -> SAV -> Stickleback ;
  DO_Temp -> SAV -> Sculpin ;
  Substrate -> Sculpin -> GOBY;
  Wind -> Zone;
  
  
  Year [exposure] ;
  DO_Temp [exposure] ;
  Zone [exposure] ;
  Stickleback [exposure] ;
  Sculpin [exposure] ;
  Microsporidia [exposure] ;
  Substrate [exposure] ;
  SAV [exposure] ;
  Breach [exposure] ;
  Rain [exposure] ;
  Wind [exposure] ;
}")

## add in zone?
par(mfrow = c(1,1))

plot(DAG_GOBY)
impliedConditionalIndependencies(DAG_GOBY)


# to pretty up the plot 
# not done yet
coordinates(DAG_GOBY) <- list(x=c(Year=1,DO_Temp=3,
                                  Stickleback=3,Sculpin=4,
                                  Rain=1,Microsporidia=1, 
                                  Substrate=3, SAV=5, 
                                  Breach=3, Wind=5, 
                                  Zone= 5, GOBY=1),
                              y=c(Year=3,DO_Temp=-3,
                                  Stickleback=-1,Sculpin=1,
                                  Rain=-5,Microsporidia=-3, 
                                  Substrate=3, SAV=-1, 
                                  Breach=-5, Wind=-5, 
                                  Zone=-3, GOBY=0))
plot(DAG_GOBY)


#####----------------------------------
#try with ggdag
library(ggdag)
goby_dag <- dagify(
  GOBY ~ Year,
  GOBY ~ O2,
  GOBY ~  Stickleback,
  GOBY ~ Sculpin,
  GOBY ~ Microsporidia,
  GOBY ~ Substrate,
  GOBY ~ SAV,
  Breach ~ Rain,
  Breach ~ Wind,
  O2 ~ Breach,
  #O2 ~ Rain,
  O2 ~ Wind,
  Sculpin ~ O2,
  Stickleback ~ O2,
  Stickleback ~ SAV,
  Sculpin ~ SAV,
  Sculpin ~ Substrate,
  
  exposure = "SAV",#"Breach",
  outcome = "GOBY"
  #labels = c(
  #  "coffee" = "Coffee", 
  #  "cancer" = "Lung Cancer", 
  #  "smoking" = "Smoking", 
  #  "addictive" = "Addictive \nBehavior"
  #)
)
ggdag(goby_dag, text_col = "white", 
      stylized = TRUE) + 
  theme_dag_blank()
ggdag_paths(goby_dag, 
            #adjust_for = c("Breach"),
            text_col = "black")

ggdag_adjustment_set(goby_dag, text_col = "black")

library(dagitty)

adjustmentSets(goby_dag)