#############################################
#  Binomial logistic regression of seed set #
#############################################

.libPaths(c("C:\\rlib", .libPaths()))
setwd("C:\\Users\\461776\\Dropbox\\PhD Hull\\Work\\Data and analysis\\Chapter 3\\Chapter-3")
dframe1 <- read.csv("Data\\SeedSetBinom.csv")
# dframe1 <- read.csv(file.choose())

names(dframe1)


dframe1$fPlot <- factor(dframe1$Plot) # treat "Plot", "Round", "PlantNo" as factors
dframe1$fRound <- factor(dframe1$Round)
dframe1$fPlantNo <- factor(dframe1$PlantNo)
dframe1$fDistance <- factor(dframe1$Distance) # make "Distance" available as a factor if reqd

dframe1$Regime <- relevel(dframe1$Regime,"Control") # relevel variables related to control level
dframe1$Pollinators <- relevel(dframe1$Pollinators,"Control")
dframe1$LitUnlit <- relevel(dframe1$LitUnlit, "Unlit")

summary(dframe1)

with(dframe1,table(Plot,Round)) # plants per plot and round - all plants included as they should be

library(lme4)  # loading up the libraries
library(car)
source("CheckResidsFunction.R") # a function for plotting glmer residuals (deviance and Pearson)

### plots

#use install.packages("ggplot2") if necessary
library(ggplot2)
source("MultiplotFunction.R") # function for panel plots in ggplot2 - see http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/


## glm for covariates

model1 <- glm(cbind(Successes, Failures) ~ Light + Regime + Pollinators + fDistance,
              family = binomial
              (link = "logit"),
              data = dframe1)

summary(model1)
drop1(model1, test = "Chi")

# fine, but does pseudoreplication of plot, round and plant have an impact (randeffs)?:

## glmm

model2 <- glmer(cbind(Successes, Failures) ~ Light + Regime + Pollinators + Distance
                +(1|fPlot) + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1)

summary(model2)
drop1(model2, test = "Chi")


# how about analysing each confounded variable within subset of data?

# so we could next examine effect of Regime, GIVEN that for a regime to exist, there must be lighting
# therefore just within AllNight + Midnight Regime treatments. Dropping the Control will create a nice fully-crossed dataset for this model...

dframe1a <- subset(dframe1,Regime!="Control") # Midnight + AllNight
summary(dframe1a)

model5 <- glmer(cbind(Successes, Failures) ~ Light + Regime + Pollinators + Distance
                + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1a)

summary(model5)
drop1(model5, test = "Chi")

chkres(model5) # fine

# great, this worked nicely and allows us to say (1) there is no sig diff between HPS and LED lights;
# and (2) there is a significantly higher rate of pollination under full-night lighting than part-night (bit weird!)

# but this doesn't answer any questions of lighting vs unlit... so what now?

# could pool lighting treatments as 'lit':
# visual inspection of full dataset suggests that this might be problematic; there is a strong risk of Type I error for Light. (But let's try and see!)


model6 <- glmer(cbind(Successes, Failures) ~ LitUnlit + Pollinators + Distance
                + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1)

summary(model6)
drop1(model6, test = "Chi") # effect of *being lit*, all else being equal, is non-sig. So no Type I error.

chkres(model6)

# just out of curiosity, what happens when you pretend regime doesn't exist...

model7 <- glmer(cbind(Successes, Failures) ~ Light + Pollinators + Distance
                + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1)

summary(model7)
drop1(model7, test = "Chi")  # EVEN with all that variation attributable to Regime left to slot into the Light variable, Light is still non-significant.

# but could that be because midnight=control masking an effect? Let's try and find a way to remove the midnight data from the model.

dframe1b <- subset(dframe1,Regime!="Midnight") # Control + AllNight
summary(dframe1b)

model8 <- glmer(cbind(Successes, Failures) ~ Light + Regime + Pollinators + Distance
                + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1)

summary(model8)
drop1(model8, test = "Chi") 

## so this is still confounded. An alternative is to analyse the two lighting types separately, with Regime coded as All, Mid, and None (=Control) in each
# first analyse HPS as the major form of lighting
library(plyr) # required for 'revalue' function

dframe1c <- subset(dframe1,Light!="LED") # HPS + CON in dataset
dframe1c$Regime <- revalue(dframe1c$Regime, c("Control"="None"))
dframe1c$Light  <- revalue(dframe1c$Light, c("CON"="HPS"))

summary(dframe1c)


#glmm
model9 <- glmer(cbind(Successes, Failures) ~ Regime + Pollinators + Distance
                + (1|fPlantNo) + (1|fRound) + (1|fPlot),
                family = binomial (link = "logit"),
                data = dframe1c)

summary(model9)
drop1(model9, test = "Chi") 


chkres(model9) # fine


# what about interaction?
model9a <- glmer(cbind(Successes, Failures) ~ Regime:Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound) + (1|fPlot),
                 family = binomial (link = "logit"),
                 data = dframe1c)

summary(model9a)
drop1(model9a, test = "Chi")

# looking solely within HPS lighting, there is a significant effect of regime on pollination (driven by AllNight);
# how about with the new challenger, LED?

dframe1d <- subset(dframe1,Light!="HPS") # LED + CON in dataset
dframe1d$Regime <- revalue(dframe1d$Regime, c("Control"="None"))
dframe1d$Light  <- revalue(dframe1d$Light, c("CON"="LED"))

summary(dframe1d)

#glmm
model10 <- glmer(cbind(Successes, Failures) ~ Regime + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound) + (1|fPlot),
                 family = binomial (link = "logit"),
                 data = dframe1d)

summary(model10)
drop1(model10, test = "Chi") 

chkres(model10) # residuals fine

# so, with LED also there is a significant effect of regime on pollination (driven by AllNight)

# interaction
model10a <- glmer(cbind(Successes, Failures) ~ Regime:Pollinators + Distance
                  + (1|fPlantNo) + (1|fRound) + (1|fPlot),
                  family = binomial (link = "logit"),
                  data = dframe1d)

summary(model10a)
drop1(model10a, test = "Chi")

## how about a comparison between LED and HPS? Let's try the same tactic and break up AllNight & Midnight

dframe1e <- subset(dframe1,Regime!="Midnight") # ALL + CON in dataset
dframe1e$Regime <- revalue(dframe1e$Regime, c("Control"="AllNight"))

summary(dframe1e)

#glmm
model11 <- glmer(cbind(Successes, Failures) ~ Light + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1e)

summary(model11)
drop1(model11, test = "Chi") 

chkres(model11)

# no effect (marginally significant effect) of light within AllNight regime;
# how about new challenger Midnight (= part-night)

dframe1f <- subset(dframe1,Regime!="AllNight") # MID + CON in dataset
dframe1f$Regime <- revalue(dframe1f$Regime, c("Control"="Midnight"))

summary(dframe1f)

#glmm
model12 <- glmer(cbind(Successes, Failures) ~ Light + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1f)

summary(model12)
drop1(model12, test = "Chi") 

chkres(model12)

# no significant effect of light within Midnight


### SO. Probability of seed set is affected by:
#(a) pollinator regime - both diurnal and nocturnal contribute (noct slightly more), with some redundancy.
#(b) lighting regime - pollination significantly boosted under full night lighting, but no difference between unlit and part-night lighting
# no difference between HPS and LED
# no difference at different distances from 0-20m

# although this is a positive effect for the plants under full-night streetlights, it is clear evidence of a disruption of overall pollination services
# raises the question of what happens on a longer transect away from the light - at what distance does effect become negative; at what distance is there no effect?

