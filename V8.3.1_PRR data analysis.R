# AUTOMATED DATA ANALYSIS FOR THE Plasmodium falciparum IN VITRO PARASITE REDUCTION RATIO ASSAY V2
# Version: 8.3.1 (21-09-2022)

# Copyright 2022 PRR assay V2 core team

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.


rm(list=ls()) ## Delete the environment
#install.packages(c("rstudioapi", "BBmisc", "plyr", "dplyr", "readxl", "ggplot2", "openxlsx")) ## You only need to install these once

### Load all packages that you will need from the library
library(rstudioapi) ## used to retrieve the path where script and data are stored
library(BBmisc)     ## used to allow for errors in compilation of the model 
library(plyr)       ## used for function 'ddply'
library(dplyr)      ## used for function 'bind_rows'
library(ggplot2)    ## used for plotting
library(readxl)     ## used to read an excel file
library(openxlsx)   ## used to create an excel file

### Specify path to data and load data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) ## Set directory to the folder containing this R script

measurements <- read_excel("data/input_example.xlsx") ## Load the excel file with the data; enter the correct file name by replacing 'input_example'
measurements$new = paste0(measurements$drug, "_", measurements$campaign) ## create a new column that creates a unique identifier which combines drug name and campaign to separate them


#### CALCULATE THE GROWTH RATE ----
controls <- measurements[measurements$drug == "control",] ## Controls only

controllist      = list()
controllistsim   = list()
controllistdata  = list()
controllistmean  = list()
controlplot      = list()
controlcounter = 1

for (campaignk in unique(controls$campaign)) {
  
  subcamp = subset(controls, controls$campaign == campaignk)
  subcamp$logP = subcamp$logP * log(10) ## growth rate is usually processed at natural logarithm scale in mathematical modeling, so we need to transform data from log10 to ln scale
  modelcamp = nls(logP ~ kgrowth * exposure + log(100000), data = subcamp, start = list(kgrowth=0.01), lower = list(kgrowth=0), upper = list(kgrowth=0.2), algorithm = "port")
  summarycamp = summary(modelcamp)
  kgrowth   = summarycamp$coefficients[1] ## growth rate of untreated control
  kgrowth5  = summarycamp$coefficients[1]  - 1.96 * summarycamp$coefficients[2]
  kgrowth95 = summarycamp$coefficients[1]  + 1.96 * summarycamp$coefficients[2]
  controlframe = data.frame("campaign" = campaignk, "kgrowth" = kgrowth, "text" = paste0(round(kgrowth,3), " [", round(kgrowth5, 3), " - ", round(kgrowth95, 3), "]"))
  
  sim = data.frame("x" = seq(0,48,1))           ## a data frame with column 'x' which is exposure time [h] with a resolution of 1 h
  sim$y = kgrowth * sim$x + log(100000)         ## use the same equation in the model to simulate the fit by applying it on 'x' 
  sim$y5 =  kgrowth5 * sim$x + log(100000)      ## simulate the bottom 5th percentile
  sim$y95 =  kgrowth95 * sim$x + log(100000)    ## simulate the top 95th percentile
  sim$campaign = campaignk                      ## add the drug_campaign name to the data frame
  
  means = ddply(.data = subcamp, .variables = c("exposure", "campaign"), summarise, "mean" = mean(logP, na.rm=TRUE)) ## Calculate the mean logP values
  
  newlab <- bquote(" " ~ Log[10] ~ "PMR" == .(round(kgrowth * 48 /log(10), 2)))
  controlplot = ggplot()+
    geom_point(data=subcamp, aes(x=exposure, y=logP), size=2,  color="black", alpha = 0.3) +      ## black points, all the observations
    geom_point(data=means, aes(x=exposure, y=mean), shape=21, size=4, fill="yellow", alpha=0.7) + ## yellow points, the means for each time point
    geom_line(data=sim, aes(x=x, y=y), color="cornflowerblue") +                                  ## blue line, the simulation
    geom_ribbon(data=sim, aes(x=x, ymin=y5, ymax=y95), fill = "cornflowerblue", alpha=0.2) +      ## ribbon, the confidence interval
    scale_x_continuous(breaks=c(0,24,48), labels=c(0,24,48), limits=c(-2,50)) +
    scale_y_continuous(breaks=seq(10, 18, 2), labels=seq(10, 18, 2), limits=c(8.5,19)) +
    xlab(paste("Exposure time [h]")) +
    ylab("Ln(viable parasites +1)") +
    annotate("text", x=0, y=16.5, label=deparse(newlab), parse=TRUE, size=3.5, hjust=0, check_overlap=T) +
    geom_text(data=controlframe, aes(x=0, y=17.5, label=paste0("  Growth rate = ", text)), size=3.5, show.legend=FALSE, hjust=0) +   ## add the parameters to the plot
    ggtitle("Growth control") +    
    theme_bw()
  
  ggsave(controlplot, filename = paste0(subcamp$campaign[1], "_growthcontrol.png"), path = "figures/", width=6, height=4) ## Create a .png file, change file name by replacing "controlplot" in "_controlplot.png"
  
  controllist[[controlcounter]]      = controlframe
  controllistdata[[controlcounter]]  = subcamp
  controllistsim[[controlcounter]]   = sim
  controllistmean[[controlcounter]]  = means
  controlcounter = controlcounter + 1
}

control = bind_rows(controllist)
data    = bind_rows(controllistdata)
sim     = bind_rows(controllistsim)
means   = bind_rows(controllistmean)


#### CALCULATE THE PHARMACODYNAMIC PARAMETERS FOR EACH COMPOUND ----
measurements <- measurements[!measurements$drug == "control",] ## Compounds only

listFLAGS      = list() ## create a new list which will contain abnormal data points, i.e. which deviate too much from fitted curve (i.e. the red circles on the plot)
listparams     = list() ## create a new list which will contain the parameters needed for plotting
listOBS        = list() ## create a new list which will contain the observations for each unique drug_campaign in a unique element (i.e. the grey/black points on the plot)
listmeans      = list() ## create a new list which will contain the means for each time point for the purposes of plotting (i.e. the yellow circles on the plot)
listparameters = list() ## create a new list which will contain the final parameters for each drug_campaign in a new element (i.e. PRR / Lag time / PCT etc.)
listsim        = list() ## create a new list which will contain the final simulation  for each drug_campaign in a new element (i.e. the blue line on the plot)
drugcounter = 1         ## create a counter that increases by 1 in every loop to populate each list


for (drugk in unique(measurements$new)) { ## cycling through each drug_campaign to calculate its parameters
  
  sub = subset(measurements, measurements$new == drugk)                                  ## create a new data frame that only contains one drug_campaign which will be defined in each loop
  drugid = sub$drug[1]
  means = ddply(.data = sub, .variables = "exposure", summarise, "mean" = mean((logP)))  ## calculate the means for each time point (i.e. the yellow circles on the plot)
  means$sd = tapply(sub$logP, sub$exposure,"sd")                                         ## calculate the standard deviation for each time point
  means$new = drugk                                                                      ## specify the unique drug identifier
  means$campaign = unique(sub$campaign)                                                  ## specify the campaign
  means$drug = drugid                                                                    ## specify the drug name
  
  data = data.frame("x" = sub$exposure, "y" = sub$logP)                                  ## create a data frame that contains time (column x) and log parasites (column y) -- made for simplicity and cleaning up
  data$new = drugk                                                                       ## add the unique drug identifier
  data$campaign = unique(sub$campaign)                                                   ## add the campaign
  data$drug = drugid                                                                     ## add the drug name
  
  meanSD = 0.5064584                                                                     ## average standard deviation (SD); pre-defined on the basis of > 40 compounds
  
  ### Tail phase
  if((subset(means$mean, means$exposure == 96)-subset(means$mean, means$exposure == 120)) > meanSD) {
    minima <- 0 
  } else { minima <- min(means$mean)} ## the mean of the minimum observed log parasites 
  
  
  ### Lag phase
  lags1 = c(0.1, 6, 12,  18, 23.9, 30, 36, 42, 47.9, 54, 60, 66, 71.9)                  ## create a list of dummy lag times use by the model
  lags = data.frame("LAG" = lags1)                                                      ## create data frame containing these lag times
  lags$sigma = Inf                                                                      ## create a column for the sigma values (residual standard deviation), the default is set to infinity
  counter = 1                                                                           ## set a counter for the coming loop that will cycle through all the lag times 
  
  for (lag in lags$LAG) { ## cycle through the dummy lag times and fit the model by fixing the lag time to each of them
    
    ## If the nls returns an error (i.e. the fit is very, very bad), sigma is set to infinity for that specific lag phase and continue
    if(is.error(try(nls(formula =  y ~ ifelse(x<=lag, 5, pmax(-a * (x-lag) + 5, minima)), start=list(a=0.05), 
                        data=data, algorithm = "port", lower=list(a=0), upper=list(a=1))))){ 
      lags$sigma[counter] = Inf  
    } else {
      
      ## If the nls does not return an error, we can safely run the model:
      model = nls(formula =  y ~ ifelse(x<=lag,                          ## 1) Lag phase: 
                                        5,                               ##    the model assumes that logP = 5 as long as x < lag (one of the pre-defined dummy lag times)           
                                        pmax(-a * (x-lag) + 5, minima)), ## 2) Linear phase (-a * (x-lag) + 5): 
                                                                         ##    when x >= lag, the model is fitted to the data. The curve changes with slope 'a' until logP reaches the pre-defined minimum 'minima'
                                                                         ## 3) Tail phase:
                                                                         ##    when logP <= minima, then logP = minima 
                  start=list(a=0.05), ## initial guess for slope 'a', the rate of parasite decline
                  data=data,          ## the data being analyzed
                  algorithm = "port", ## algorithm port is preferred as it allows setting lower and upper parameters
                  lower=list(a=0),    ## parasite decline not allowed to be less than zero
                  upper=list(a=1))    ## parasite decline not allowed to be higher than 1 (i.e. PRR = 48)
      
    if(is.error(try(summary(model)))) {   ## if the model summary returns an error (i.e. the fit is very, very bad), sigma is set to infinity
        lags$sigma[counter] = Inf         
      } else {                            ## if the model summary does does not return an error, we can retrieve the sigma value from it
        summary = summary(model)            
        lags$sigma[counter] = summary$sigma 
      }}
    counter=counter+1                     ## add 1 to the counter for each assessed lag time so we carry on populating the lags data frame
  }
  
  
  ### Dominant rules for lag phase determination 
  
  ## 1) if one of the "real" lag times has a sigma <= 1% higher than the lowest sigma, prefer this lag time (or the one with the lowest sigma if more than one real lag times are good)
  lags$diff = NA                                  
  for (i in 1:(nrow(lags))) {                     
    lags$diff[i] = lags$sigma[i]/min(lags$sigma)  
  }
  lags$diffabs = abs(1-lags$diff) ## an absolute measure of the difference between each lag. for example the previous step can result in 0.99 or 1.01 , both are very similar so take the absolute difference between that and 1 and this will tell us the difference whether it's more or less
  realms <- c(0.1, 23.9, 47.9, 71.9) ## a vector containing only the real measurement times
  truelags <- subset(lags, lags$LAG %in% realms & lags$diffabs <=0.01) ## check if one of the real lag times produces a sigma which is <= 1% different from lag time with minimum sigma value
  lags0 <- lags
  
  if(nrow(truelags) > 0) {  
    lags = truelags
  } 
  
  if(!all(is.infinite(lags$sigma))) { ## if not all sigma are infinity (i.e. we have at least one run that converged)
    lag = lags$LAG[which.min(lags$sigma)] ## pick the lag time that returns the least sigma value
    
  ## 2) if at time x, logP + 0.5* average SD >= 5, lag phase >= x
    if ((means$mean[means$exposure == 24] + 0.5*meanSD) >= 5 & lag < 23.9) { # for x=24h
      lags3 <- subset(lags0, lags0$LAG >= 23.9) 
      lag <- lags3$LAG[which.min(lags3$sigma)]}
    
    if ((means$mean[means$exposure == 48] + 0.5*meanSD) >= 5 & lag < 47.9) { # for x=48h
      lags3 <- subset(lags0, lags0$LAG >= 47.9)  
      lag <- lags3$LAG[which.min(lags3$sigma)]}
    
    if ((means$mean[means$exposure == 72] + 0.5*meanSD) >= 5 & lag < 71.9) { # for x=72h
      lags3 <- subset(lags0, lags0$LAG >= 71.9)  
      lag <- lags3$LAG[which.min(lags3$sigma)]}  
    
  ## 3) if at time x, with x>24h, logP + 1* average SD >= 5, lag phase >= x-12
    if ((means$mean[means$exposure == 48] + meanSD) >= 5 & lag < 23.9) { # for x=48h
      lags4 <- subset(lags0, lags0$LAG >= 36)
      lag <- lags4$LAG[which.min(lags4$sigma)]}
    
    if ((means$mean[means$exposure == 72] + meanSD) >= 5 & lag < 47.9) { # for x=72h
      lags4 <- subset(lags0, lags0$LAG >= 60)
      lag <- lags4$LAG[which.min(lags4$sigma)]}
    
  ## 4) if logP drop between 0h and 24h >=1.5, lag phase < 12
          # if logP drop between 0h and 24h >=2,   lag phase = 0 
          # if logP drop between time x and x+24h >=1.5 and lag >x, lag phase = x+6 (for x=24 or x=48)
          # if logP drop between time x and x+24h >=2   and lag >x, lag phase = x   (for x=24 or x=48)
    
    if ((means$mean[means$exposure == 0] - means$mean[means$exposure == 24]) >= 1.5) { # for 0-24h window
      lags2 <- subset(lags0, lags0$LAG < 12) 
      lag <- lags2$LAG[which.min(lags2$sigma)]}
    
    if ((means$mean[means$exposure == 0] - means$mean[means$exposure == 24]) >= 2) {
      lag = 0.1}
    
    if ((means$mean[means$exposure == 24] - means$mean[means$exposure == 48]) >= 1.5 & lag > 23.9) { # for 24-48h window
      lag = 30} 
    
    if ((means$mean[means$exposure == 24] - means$mean[means$exposure == 48]) >= 2 & lag > 23.9) { 
      lag = 23.9}
    
    if ((means$mean[means$exposure == 48] - means$mean[means$exposure == 72]) >= 1.5 & lag > 47.9) { # for 48-72h window
      lag = 54}
    
    if ((means$mean[means$exposure == 48] - means$mean[means$exposure == 72]) >= 2 & lag > 47.9) { 
      lag = 47.9}
    
    
    ### Linear phase 
    ### Run the nls function with the pre-defined tail- and lag phase
    
    model = nls(formula =  y ~ ifelse(x<=lag, 5, pmax(-a * (x-lag) + 5, minima)), start=list(a=0.05), data=data, algorithm = "port", lower=list(a=0), upper=list(a=1)) 
    summary = summary(model) ## define the summary of the model
    
    a = summary$coefficients[1]    ## define the estimated parameter 'a', the rate of parasite decline
    aSE = summary$coefficients[2]  ## the standard error of 'a'
    a5 = a - (1.96* aSE)           ## the lower/5th percentile as it is defined mathematically
    a95 = a + (1.96* aSE)          ## upper/95th percentile as it is defined mathematically
    
    ## Calculate PRR within 48h (+ 95% CI)
    PRR = a*48 
    PRR5 = a5*48
    PRR95 = a95*48
    
    ## Calculate PCT (+ 95% CI)
    PCT99.9_a = 3/a + lag    
    PCT99.9_5 = 3/a95 + lag  
    PCT99.9_95 = 3/a5 + lag  
    
    ## Calculate Emax (+ 95% CI)
    ## Note that Emax only equals the maximum killing effect of the drug assuming that the chosen concentration corresponds to a concentration at which killing has saturated
    Kgrowth = control$kgrowth[control$campaign == unique(sub$campaign)] ## growth rate for this specific campaign
    
    Emax = a * log(10) + Kgrowth      ## Calculate Emax by adding 'a' to the growth rate. Note: Need to transform 'a' into ln scale first, so that Emax is at ln scale, too.
    Emax5 = a5 * log(10) + Kgrowth    
    Emax95 = a95 * log(10) + Kgrowth  
    
    ## Define the P value for the fit 
    P = summary$coefficients[4] 
    
    ### Define the pharmacodynamic category of each compound (slow, intermediate, fast)
    if(round(lag, 0) == 0 && PRR >= 4){category <- "fast"}
    if(round(lag, 0) > 0  && PRR >= 4){category <- "fast with lag phase"}
    if(round(lag, 0) == 0 && 3 <= PRR && PRR < 4){category <- "intermediate"}
    if(round(lag, 0) > 0  && 3 <= PRR && PRR < 4){category <- "intermediate with lag phase"}
    if(round(lag, 0) > 48 || PRR < 3){category <- "slow"}
    
  } else { ## if all tested lag times returned an error (i.e. all sigmas are infinity) define all the parameters as 'NA' 
    a = NA; aSE = NA; a5 = NA; a95 = NA
    PRR = NA; PRR5 = NA; PRR95 = NA
    PCT99.9 = NA; PCT99.9_5 = NA; PCT99.9_95 = NA 
    Emax = NA; Emax5 = NA; Emax95 = NA; P = NA
    category = "NA"; lag = 120
  }
  
  
  ### Create a simulation based on the estimated parameters for each drug campaign
  sim = data.frame("x" = seq(0,24*5,1))                                 ## a data frame with column x which is time in hours from 0 to 120 with a  resolution of 1 hour
  sim$y   = ifelse(sim$x<=lag, 5, pmax(-a * (sim$x-lag) + 5, minima))   ## use the same equation in the model to simulate the fit by applying it on the times ('a' is now already estimated and the estimated value is used)
  sim$y5  = ifelse(sim$x<=lag, 5, pmax(-a5 * (sim$x-lag) + 5, minima))  ## simulate the bottom 5th percentile
  sim$y95 = ifelse(sim$x<=lag, 5, pmax(-a95 * (sim$x-lag) + 5, minima)) ## simulate the top 95th percentile
  sim$new = drugk                                                       ## add the unique identifier to the data frame
  sim$drug = drugid                                                     ## add the drug name to the data frame
  sim$campaign = unique(sub$campaign)                                   ## add the campaign ID to the data frame
  
  PCT99.9 = ifelse(min(sim$y)>2, NA, PCT99.9_a) ## if logP never falls below 2, the PCT cannot be determined and is set to 'NA'
  PRRclinical <- round(as.numeric(subset(means, exposure == 0, mean) - subset(means, exposure == 48, mean)), 1) ## calculate PRR clinical = mean logP difference between 0 and 48h
  
  ### Flag suspicious data points in the graph (e.g. with large standard deviation, with large deviation from simulation)
  means$sim = subset(sim$y, sim$x %in% c(0, 24, 48, 72, 96, 120))
  highlight = means[means$sd > 2*meanSD | abs(means$mean-means$sim) > meanSD, ] 
  highlight = unique(highlight)
    comment = as.character(); comment = ifelse(nrow(highlight) == 0,  "NA", "potential abnormality detected") ## warning that will appear in the results table
  
  
  ### Create a table that has the parameters that were calculated above
  parameters = data.frame("PRR" = paste0(round(PRR,1), " [", round(PRR5, 1), " - ", round(PRR95, 1), "]"),
                          "PCT99.9" = ifelse(is.na(PCT99.9), ">120 h",   paste0(round(PCT99.9,1), " [", round(PCT99.9_5, 1), " - ", round(PCT99.9_95, 1), "] h")),
                          "LagTime" = as.numeric(round(lag,0)),
                          "Emax" = paste0(round(Emax,2), " [", round(Emax5, 2), " - ", round(Emax95, 2), "] /h"),
                          "Kgrowth" = subset(control$text, control$campaign == unique(sub$campaign)), 
                          "new" = drugk, 
                          "campaign" = unique(sub$campaign)) 
  
  ### Create a table that has the parameters that were calculated above for the results table
  params = data.frame("Drug ID" = unique(sub$drug), 
                      "campaign" = unique(sub$campaign), 
                      "hill slope" = round(as.numeric(unique(sub$hill)), 3),
                      "IC50 [uM]" = round(as.numeric(unique(sub$IC50)), 4), 
                      "multiple of IC50" = unique(sub$factor), 
                      "category" = category, 
                      "lag phase [h]" = ifelse(round(lag, 0) %in% c(0, 24, 48, 72), paste0(round(lag, 0), ""), paste0(round(lag,0), " [", round_any(lag, 24, floor), " - ", round_any(lag, 24, ceiling), "]")), 
                      "99.9% PCT [h]" = ifelse(is.na(PCT99.9), ">120",   paste0(round(PCT99.9,1), " [", round(PCT99.9_5, 1), " - ", round(PCT99.9_95, 1), "]")),
                      "PRR" = paste0(round(PRR,1), " [", round(PRR5, 1), " - ", round(PRR95, 1), "]"),
                      "Emax [/h]" = paste0(round(Emax,2), " [", round(Emax5, 2), " - ", round(Emax95, 2), "]"),
                      "growth rate [/h]" = subset(control$text, control$campaign == unique(sub$campaign)),
                      "sigma" = round(min(summary$sigma), 3), 
                      "PRR clinical" = PRRclinical,
                      "comment" = comment,
                      check.names = F)
  
  ### Start populating the empty lists that were created
  listOBS[[drugcounter]]        = data 
  listmeans[[drugcounter]]      = means
  listparameters[[drugcounter]] = parameters
  listsim[[drugcounter]]        = sim
  listparams[[drugcounter]]     = params
  listFLAGS[[drugcounter]]      = highlight
  
  ### Create a plot for each individual compound of the campaign
  oneplot <- ggplot() + 
    geom_line(data=sim, aes(x=x, y=y), color="cornflowerblue") +
    geom_ribbon(data=sim, aes(x=x, ymin=y5, ymax=y95), fill = "cornflowerblue", alpha=0.2) + 
    geom_point(data=sub, aes(x=exposure, y=logP), size=2,  color="black", alpha = 0.3) + 
    geom_point(data=means, aes(x=exposure, y=mean), shape=21, size=4, fill="gold", alpha=0.7) + 
    geom_point(data = highlight, aes(x=exposure, y=mean), shape=21, size=4, fill="firebrick", alpha=0.7) +
    xlab(paste("Exposure time [h]")) +
    ylab(expression(paste(Log[10], "(viable parasites +1)"))) +
    scale_x_continuous(breaks=seq(0,120,24), limits=c(-5,125)) +
    scale_y_continuous(breaks=seq(0,6,1), limits=c(-0.25,6)) +
    ggtitle(sub$drug) +
    theme(panel.grid.minor=element_blank(), axis.title = element_text( size=11), panel.background = element_blank(), panel.grid.major = element_line(color="grey", size=0.25), axis.line=element_line(color="grey"),
          axis.ticks = element_line(colour = "grey"), strip.text = element_text(size=9), axis.text=element_text(size=9),   panel.border = element_rect(colour = "grey", fill=NA, size=1), plot.title = element_text(size = 12))
  
   ggsave(oneplot, filename = paste(campaignk, "_", drugid, ".png"), path = "figures/", width=3.2, height=3.2) ## save as .png file 
  
  drugcounter = drugcounter+1 ## end of the loop; add 1 to the drug counter so next time the lists are populated appropriately
}


## convert the big lists that were created before into big tables that will be used for plotting
theflags      = bind_rows(listFLAGS)
thedata       = bind_rows(listOBS)
themeans      = bind_rows(listmeans)
theparameters = bind_rows(listparameters)
theparams     = bind_rows(listparams)
thesim        = bind_rows(listsim)


#### CREATE A TABLE CONTAINING THE PHARMACODYNAMIC PARAMETERS AND A SUMMARY PLOT FOR EACH CAMPAIGN ----

for (campaignK in unique(measurements$campaign)) { ## Cycle through each campaign separately
  
  campmeasurements = subset(measurements, measurements$campaign == campaignK)
  campdata = subset(thedata, thedata$campaign == campaignK)
  campmeans = subset(themeans, themeans$campaign == campaignK)
  campflags = subset(theflags, theflags$campaign == campaignK)
  campsim = subset(thesim, thesim$campaign == campaignK)
  campparameters = subset(theparameters, theparameters$campaign == campaignK)
  campparams = subset(theparams, theparams$campaign == campaignK)
  
  ### Crate a summary plot for whole campaign (.png)
  theplot = ggplot()+
    geom_point(data=campdata,  aes(x=x, y=y), color="black", size=2,  alpha=0.3) + ## Black points, all the observations
    geom_point(data=campmeans, aes(x=exposure, y=mean),shape=21, size=4, fill="yellow", alpha=0.7) + ## Yellow points, the means for each time point
    geom_point(data=campflags, aes(x=exposure, y=mean),shape=21, size=4, fill="firebrick", alpha=0.7) +
    geom_line(data=campsim,    aes(x=x, y=y), color="black") + ## Black line, the simulation
    geom_ribbon(data=campsim,  aes(x=x, ymin=y5, ymax=y95), fill="cornflowerblue", alpha=0.2) + ## Blue ribbon, 95% confidence interval
    xlab(paste("Exposure time [h]")) +
    ylab(expression(paste(Log[10], "(viable parasites +1)"))) +
    scale_x_continuous(breaks=seq(0,120,24), limits=c(-5,125)) +
    scale_y_continuous(breaks=seq(0,6,1),    limits=c(-0.25,6)) +
    geom_text(data=campparameters, aes(x=50, y=6,   label=paste0("PRR = ", PRR)),         size=3, show.legend=FALSE, hjust=0) + ## Add the parameters to the plot
    geom_text(data=campparameters, aes(x=50, y=5.5, label=paste0("PCT99.9 = ", PCT99.9)), size=3, show.legend=FALSE, hjust=0) +
    geom_text(data=campparameters, aes(x=50, y=5.0, label=paste0("Emax = ", Emax)),       size=3, show.legend=FALSE, hjust=0) +
    geom_text(data=campparameters, aes(x=50, y=4.5, label=paste0("Lag Time = ",LagTime)), size=3, show.legend=FALSE, hjust=0) +
    geom_text(data=campparameters, aes(x=50, y=4.0, label=paste0("Kgrowth = ", Kgrowth)), size=3, show.legend=FALSE, hjust=0) +
    facet_wrap(~new, scales = "free")+ ##Facet by the drug_campaign 
    theme(panel.grid.minor=element_blank(), axis.title=element_text(face="bold", size=12), panel.background=element_blank(), panel.grid.major=element_line(color="grey", size=0.25), 
          axis.line=element_line(color="black"), strip.text=element_text(size=10, face="bold"), axis.text=element_text(size=10))
  
  ggsave(theplot, filename=paste(campaignK, "_summary.png"), path="figures/", width=11, height=6) ## Save all plots in one .png file (for 1 campaign / 6 compounds)
  
  
  ### Create a result table containing the pharmacodynamic parameters of the whole campaign (.xlsx)
  wb2 = createWorkbook()
  addWorksheet(wb2, "results")
  writeData(wb2, "results", campparams, startRow = 1, startCol = 1); setColWidths(wb2, sheet = "results", cols = 1:ncol(campparams), widths = "auto") 
  saveWorkbook(wb2, paste0("results/", campaignK, "_results.xlsx")) ## Change Excel file name here
}








