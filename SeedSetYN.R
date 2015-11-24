#################################
### Silene latifolia seed set ###
#################################

.libPaths(c("C:\\rlib", .libPaths()))
setwd("C:\\Users\\461776\\Dropbox\\PhD Hull\\Work\\Data and analysis\\Chapter 3\\Chapter-3")
dframe1 <- read.csv("Data\\SeedSetYN.csv")
# dframe1 <- read.csv(file.choose())

names(dframe1)

dframe1$fPlot <- factor(dframe1$Plot) # treat "Plot", "Round", "PlantNo" as factors
dframe1$fRound <- factor(dframe1$Round)
dframe1$fPlantNo <- factor(dframe1$PlantNo)
dframe1$fDistance <- factor(dframe1$Distance) # make "Distance" available as a factor if reqd

dframe1$Regime <- relevel(dframe1$Regime,"Control") # relevel variables related to control level
dframe1$Pollinators <- relevel(dframe1$Pollinators,"Control")

summary(dframe1)

with(dframe1,table(Plot,Round)) # seedheads per plot and round - fairly even spread, no obvious patterns, good

#use install.packages("lme4"), install.packages("car") if necessary
library(lme4)  # loading up the libraries
library(car)

### plots

#use install.packages("ggplot2") if necessary
library(ggplot2)
source("MultiplotFunction.R") # function for panel plots in ggplot2 - see http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

# figure Light
p1 <- ggplot(dframe1,aes(x=Light,y=SeedSetYN))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

# figure Regime
p2 <- ggplot(dframe1,aes(x=Regime,y=SeedSetYN))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

# figure Treatment - Light and Regime combined
p2a <- ggplot(dframe1,aes(x=Treatment,y=SeedSetYN))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

# figure Pollinators
p3 <- ggplot(dframe1,aes(x=Pollinators,y=SeedSetYN))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

# figure Distance
p4 <- ggplot(dframe1,aes(x=Distance,y=SeedSetYN))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

multiplot(p1, p2, p3, p4, cols=2)
p2a


### analysis

## glmm

model1 <- glmer(SeedSetYN ~ Treatment + Pollinators + Distance
                + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1)

summary(model1)
drop1(model1, test = "Chi")
Anova(model1, test = "Chi")

relgrad <- with(model4@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))

# better convergence (now less than 0.001 so should be fine), but this model says treatment only marginally significant with LRT.
# (Incidentally, significant with Type II Wald, but this is a worse test, and no good reason to use it here)
# Given confounding of light & regime, is LRT appropriate? Is there another way to construct the model?


