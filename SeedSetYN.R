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

d1 <- qplot(Light, SeedSetYN, data=dframe1)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

# figure Regime
p2 <- ggplot(dframe1,aes(x=Regime,y=SeedSetYN))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d2 <- qplot(Regime, SeedSetYN, data=dframe1)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

# figure Treatment - Light and Regime combined
p3 <- ggplot(dframe1,aes(x=Treatment,y=SeedSetYN))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d3 <- qplot(Treatment, SeedSetYN, data=dframe1)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

# figure Pollinators
p4 <- ggplot(dframe1,aes(x=Pollinators,y=SeedSetYN))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d4 <- qplot(Pollinators, SeedSetYN, data=dframe1)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

# figure Distance
p5 <- ggplot(dframe1,aes(x=Distance,y=SeedSetYN))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d5 <- qplot(Distance, SeedSetYN, data=dframe1)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

multiplot(p1, d1, p2, d2, p3, d3, p4, d4, p5, d5, cols=5)


### analysis

## nested anova

ano1 <- aov(SeedSetYN ~ Light*Regime, dframe1)
summary(ano1)

ano2 <- aov(SeedSetYN ~ Regime/Light, dframe1)
summary(ano2)

## glm

model1 <- glm(SeedSetYN ~ Light + Regime + Pollinators + fDistance,
              family = binomial
              (link = "logit"),
              data = dframe1)

summary(model1)
drop1(model1, test = "Chi")

# fine, but does pseudoreplication of plot, round and plant have an impact?:

## glmm

model2 <- glmer(SeedSetYN ~ Light + Regime + Pollinators + Distance
                +(1|fPlot) + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1)

summary(model2)
drop1(model2, test = "Chi")

# nice, but model is rank deficient due to confounding of light and regime. Try combining them?:

model3 <- glmer(SeedSetYN ~ Treatment + Pollinators + Distance
                +(1|fPlot) + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1)

summary(model3)
drop1(model3, test = "Chi")

# model failed to converge - recheck convergence with:

relgrad <- with(model3@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))

# convergence =~ 0.001 - not a huge problem but inspect rand effs - fPlot has very low variance so try removing

model4 <- glmer(SeedSetYN ~ Treatment + Pollinators + Distance
                + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1)

summary(model4)
drop1(model4, test = "Chi")
Anova(model4, test = "Chi")

relgrad <- with(model4@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))



# better convergence (now less than 0.001 so should be fine), but this model says treatment only marginally significant with LRT.
# (Incidentally, significant with Type II Wald, but this is a worse test, and no good reason to use it here)
# Given confounding of light & regime, LRT appropriate? Is there another way to construct the model?


# how about analysing each confounded variable within subset of data?

# so we could next examine effect of Regime, GIVEN that for a regime to exist, there must be lighting
# therefore just within AllNight + Midnight Regime treatments. Dropping the Control will create a nice fully-crossed dataset for this model...

dframe1a <- subset(dframe1,Regime!="Control") # Midnight + AllNight
summary(dframe1a)

model5 <- glmer(SeedSetYN ~ Light + Regime + Pollinators + Distance
                + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1a)

summary(model5)
drop1(model5, test = "Chi")

relgrad <- with(model5@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) # same issues with fPlot randeff on this model so exclude again

# great, this worked nicely and allows us to say (1) there is no sig diff between HPS and LED lights;
# and (2) there is a significantly higher rate of pollination under full-night lighting than part-night

# but this doesn't answer any questions of lighting vs unlit... so what now?

# could pool lighting treatments as 'lit' for each variable, but:
# visual inspection of full dataset suggests that this might be problematic; there is a strong risk of Type I error for Light. (But let's try and see!)

dframe1$LitUnlit <- recode(dframe1$Light, "c('HPS','LED')='Lit'; else='Unlit'")

# figure LitUnlit
p6 <- ggplot(dframe1,aes(x=LitUnlit,y=SeedSetYN))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d6 <- qplot(LitUnlit, SeedSetYN, data=dframe1)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

multiplot(p6, d6, cols=2)

model6 <- glmer(SeedSetYN ~ LitUnlit + Pollinators + Distance
                + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1)

summary(model6)
drop1(model6, test = "Chi") #effect of *being lit*, all else being equal, is non-sig'ly (v slightly) positive
