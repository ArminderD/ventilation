rm(list=ls())

# Install these packages first and then load with the below code
library(boot)
library(deSolve)
library(FME)
library(ggplot2)

# Load the saved dataframe
data_df <- read.csv("./data/ClinicA.csv", stringsAsFactors = FALSE)
data_df$Room <- NULL #remove first column which only contains the room name in first cell
# format times (to give seconds)
data_df$Time <- vapply(data_df$Time, function(x) as.numeric(strptime(x, format = "%Y-%m-%dT%H:%M:%SZ")), numeric(1))
data_df$Time <- data_df$Time - data_df$Time[1]
colnames(data_df) <- c("time", "C_out", "C_in1", "C_in2", "C_in3", "I", "Adults", "Children", "Infants")
# convert ppm to litres of CO2 per litres of air (by dividing by 1x10^6)
data_df[,c("C_out", "C_in1", "C_in2", "C_in3")] <- data_df[,c("C_out", "C_in1", "C_in2", "C_in3")] / 1e6
# take mean of three C_in measurements
data_df$C_in <- (data_df$C_in1 + data_df$C_in2 + data_df$C_in3) / 3

# CO2 generation rate in litres per second, estimated from Persily et al. (change if needed)
V <- 162.791287 * 1000  # Convert cubic m to litres by multiplying by 1000 and use volume measurements (must change acrroding to room)

# CO2 generation rate in litres per second, estimated from Persily et al. (change if needed) - mean taken for both sexes
# 1.2 MET
G_adult <- 0.00377
G_child <- 0.001975
G_infant <- 0.00105

# # 1.0 METS
# G_adult <- 0.0031
# G_child <- 0.0016
# G_infant <- 0.0009
# 
# # 1.4 METS
# G_adult <- 0.0044
# G_child <- 0.0023
# G_infant <- 0.0013
# 
# # 1.6 METS
# G_adult <-0.0050
# G_child <- 0.0026
# G_infant <- 0.0014
# 
# 
# # Mean METS
# G_adult <-0.004075
# G_child <- 0.002125
# G_infant <- 0.001175


##################################
#### To fit steady state model####
##################################

# Use linear interpolation (curve fitting) for C_out (outside level of CO_2) over time, then take the average of the 
# linear interpolation over time to get the mean value of C_out:
C_out_func <- approxfun(data_df[,c("time", "C_out")], method = "linear", rule = 2) # linear interpolation
t_vec_fine <- seq(0, max(data_df$time), by = 60) # calculate C_out every 60 seconds
C_out <- C_out_func(t_vec_fine) # calculate C_out every 60 seconds
C_out <- mean(C_out) # take the mean

## Obtain initial guess for Q by linear regression for steady-state model (since Q=IG/Cin,ss-Cout):
data_df$IG_adults <- data_df$Adults * G_adult 
data_df$IG_children <- data_df$Children* G_child 
data_df$IG_infants <- data_df$Infants * G_infant 
# data_df$IG  <- (data_df$IG_adults )
data_df$IG  <- (data_df$IG_adults + data_df$IG_children + data_df$IG_infants) # the total generation rate of CO2 by the people in the room
#data_df$IG <- data_df$I * G # the total generation rate of CO2 by the people in the room
data_df$diff_C <- data_df$C_in - C_out # difference between inside and outside CO2 concentrations
# Note that diff_C for a given time is the difference between the mean inside CO2 level for the 
# three sensors at that time, and the mean outside CO2 level across all time (C-out)
# if we plot diff_C on the y-axis and IG on the x-axis we should get a straight line with slope 1/Q
# to check that the plot actually looks like a straight line

# Constrained linear regression with intercept = 0
lm1 <- lm(formula = IG ~ diff_C + 0, data = data_df) 
# Absolute ventilation Q (l/s):
Q <- unname(lm1$coefficients)
# 95% CI of Q:
Q_conf <- confint(lm1)

# If needed, Air Changes Per Hour (ACH):
ACH <- (Q*3.6)/(V/1000)
# 95% of ACH:
ACH_conf <- (Q_conf*3.6)/(V/1000)

