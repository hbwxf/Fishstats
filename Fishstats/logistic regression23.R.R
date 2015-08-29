############################################################################
老王又修改了
sdsfwewwew
ewewewew
# Eel fun, using the eel.dat data set from http://www.sagepub.com/dsur/study/default.htm
#
setwd("~/Documents/R/logistic regression")
#
# Get the eel data set
eelData <- read.delim("eel.dat", header = TRUE)
#
str(eelData)
head(eelData)
#
# reset the factor levels for the dichotomous variables "Intervention" and 
# "Cured" to the opposite of what R defaulted them to be 
eelData$Cured <- relevel(eelData$Cured, "Not Cured")
eelData$Intervention <- relevel(eelData$Intervention, "No Treatment")
str(eelData)
# 
