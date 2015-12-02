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
source("CheckResidsFunction.R") # a function for plotting glmer residuals (deviance and Pearson)

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

# regime appears important, light not, but what about covariates and randeffs...


## glm for covariates

model1 <- glm(SeedSetYN ~ Light + Regime + Pollinators + fDistance,
              family = binomial
              (link = "logit"),
              data = dframe1)

summary(model1)
drop1(model1, test = "Chi")

# fine, but does pseudoreplication of plot, round and plant have an impact (randeffs)?:

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
glht(model4, mcp(Treatment="Tukey"))

relgrad <- with(model4@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))



# better convergence (now less than 0.001 so should be fine), but this model says treatment only marginally significant with LRT.
# (Incidentally, significant with Type II Wald, but this is a worse test, and no good reason to use it here)
# Given confounding of light & regime, is there another way to construct the model?


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

chkres(model5) # ok, not amazing, but not final model so not going to worry right now!

# great, this worked nicely and allows us to say (1) there is no sig diff between HPS and LED lights;
# and (2) there is a significantly higher rate of pollination under full-night lighting than part-night (bit weird!)

# but this doesn't answer any questions of lighting vs unlit... so what now?

# could pool lighting treatments as 'lit':
# visual inspection of full dataset suggests that this might be problematic; there is a strong risk of Type I error for Light. (But let's try and see!)

dframe1$LitUnlit <- recode(dframe1$Light, "c('HPS','LED')='Lit'; else='Unlit'")
dframe1$LitUnlit <- relevel(dframe1$LitUnlit, "Unlit")

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
drop1(model6, test = "Chi") #effect of *being lit*, all else being equal, is non-sig. So no Type I error.

chkres(model6)

# just out of curiosity, what happens when you pretend regime doesn't exist...

model7 <- glmer(SeedSetYN ~ Light + Pollinators + Distance
                + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1)

summary(model7)
drop1(model7, test = "Chi")  # EVEN with all that variation attributable to Regime left to slot into the Light variable, Light is still non-significant.

# but could that be because midnight=control masking an effect? Let's try and find a way to remove the midnight data from the model.

dframe1b <- subset(dframe1,Regime!="Midnight") # Control + AllNight
summary(dframe1b)

model8 <- glmer(SeedSetYN ~ Light + Regime + Pollinators + Distance
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

# figure within HPS
p7 <- ggplot(dframe1c,aes(x=Regime,y=SeedSetYN))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d7 <- qplot(Regime, SeedSetYN, data=dframe1c)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

multiplot(p7, d7, cols=2)

#glmm
model9 <- glmer(SeedSetYN ~ Regime + Pollinators + Distance
                + (1|fPlantNo) + (1|fRound) + (1|fPlot),
                family = binomial (link = "cauchit"),
                data = dframe1c)

summary(model9)
drop1(model9, test = "Chi") 

relgrad <- with(model9@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) # model here is fine with fPlot in

chkres(model9) # kind of bimodal. worth trying another link function? Done - logit is best option (and not too bad really)
# bit of a pattern in the binned plot but generally don't need to worry with Bernoulli GLM


# looking solely within HPS lighting, there is a significant effect of regime on pollination (driven by AllNight);
# how about with the new challenger, LED?

dframe1d <- subset(dframe1,Light!="HPS") # LED + CON in dataset
dframe1d$Regime <- revalue(dframe1d$Regime, c("Control"="None"))
dframe1d$Light  <- revalue(dframe1d$Light, c("CON"="LED"))

summary(dframe1d)

# figure within HPS
p8 <- ggplot(dframe1d,aes(x=Regime,y=SeedSetYN))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d8 <- qplot(Regime, SeedSetYN, data=dframe1d)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

multiplot(p8, d8, cols=2)

#glmm
model10 <- glmer(SeedSetYN ~ Regime + Pollinators + Distance
                + (1|fPlantNo) + (1|fRound) + (1|fPlot),
                family = binomial (link = "logit"),
                data = dframe1d)

summary(model10)
drop1(model10, test = "Chi") 

relgrad <- with(model10@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) # model here is fine with fPlot in

chkres(model10) # residuals fine

# so, with LED also there is a significant effect of regime on pollination (driven by AllNight)


## how about a comparison between LED and HPS? Let's try the same tactic and break up AllNight & Midnight

dframe1e <- subset(dframe1,Regime!="Midnight") # ALL + CON in dataset
dframe1e$Regime <- revalue(dframe1e$Regime, c("Control"="AllNight"))

summary(dframe1e)

# figure within AllNight
p9 <- ggplot(dframe1e,aes(x=Light,y=SeedSetYN))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d9 <- qplot(Light, SeedSetYN, data=dframe1e)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

multiplot(p9, d9, cols=2)

#glmm
model11 <- glmer(SeedSetYN ~ Light + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1e)

summary(model11)
drop1(model11, test = "Chi") 

relgrad <- with(model11@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) # model here does not converge with fPlot in - so removed (fine thereafter)

chkres(model11)

# no effect (marginally significant effect) of light within AllNight regime;
# how about new challenger Midnight (= part-night)

dframe1f <- subset(dframe1,Regime!="AllNight") # MID + CON in dataset
dframe1f$Regime <- revalue(dframe1f$Regime, c("Control"="Midnight"))

summary(dframe1f)

# figure within Midnight
p10 <- ggplot(dframe1f,aes(x=Light,y=SeedSetYN))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d10 <- qplot(Light, SeedSetYN, data=dframe1f)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

multiplot(p10, d10, cols=2)

#glmm
model12 <- glmer(SeedSetYN ~ Light + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1f)

summary(model12)
drop1(model12, test = "Chi") 

relgrad <- with(model12@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) # model here does not converge with fPlot in - so removed (fine thereafter)

chkres(model12)

# no significant effect of light within Midnight

### SO. Probability of seed set is affected by:
#(a) pollinator regime - both diurnal and nocturnal contribute (noct slightly more), with some redundancy.
#(b) lighting regime - pollination significantly boosted under full night lighting, but no difference between unlit and part-night lighting
# no difference between HPS and LED
# no difference at different distances from 0-20m

# although this is a positive effect for the plants under full-night streetlights, it is clear evidence of a disruption of overall pollination services
# raises the question of what happens on a longer transect away from the light - at what distance does effect become negative; at what distance is there no effect?

## next up - are there any effects on quality of pollination (seed count and weight)...

### Quality of pollination ###

dframe2 <- read.csv("Data\\SeedCount.csv")
# dframe2 <- read.csv(file.choose())

names(dframe2)

dframe2$fPlot <- factor(dframe2$Plot) # treat "Plot", "Round", "PlantNo" as factors
dframe2$fRound <- factor(dframe2$Round)
dframe2$fPlantNo <- factor(dframe2$PlantNo)
dframe2$fDistance <- factor(dframe2$Distance) # make "Distance" available as a factor if reqd

dframe2$Regime <- relevel(dframe2$Regime,"Control") # relevel variables related to control level
dframe2$Pollinators <- relevel(dframe2$Pollinators,"Control")

summary(dframe2)

with(dframe2,table(Plot,Round)) # seedheads per plot and round - fairly even spread, no obvious patterns other than round 1 very poor, good

### plots

hist(dframe2$SeedWeight)
hist(dframe2$SeedCount) # check variable distributions - both appear to be Poisson (seed count possibly overdispersed - check this)

# SeedCount

# figure Light
p11 <- ggplot(dframe2,aes(x=Light,y=SeedCount))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d11 <- qplot(Light, SeedCount, data=dframe2)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

# figure Regime
p12 <- ggplot(dframe2,aes(x=Regime,y=SeedCount))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d12 <- qplot(Regime, SeedCount, data=dframe2)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

# figure Treatment - Light and Regime combined
p13 <- ggplot(dframe2,aes(x=Treatment,y=SeedCount))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d13 <- qplot(Treatment, SeedCount, data=dframe2)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

# figure Pollinators
p14 <- ggplot(dframe2,aes(x=Pollinators,y=SeedCount))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d14 <- qplot(Pollinators, SeedCount, data=dframe2)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

# figure Distance
p15 <- ggplot(dframe2,aes(x=Distance,y=SeedCount))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d15 <- qplot(Distance, SeedCount, data=dframe2)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

multiplot(p11, d11, p12, d12, p13, d13, p14, d14, p15, d15, cols=5) # no obvious trends, CIs seem to overlap in most cases. HPS vs LED possibly?


### analysis

## nested anova

ano3 <- aov(SeedCount ~ Light*Regime, dframe2)
summary(ano3)

ano4 <- aov(SeedCount ~ Regime/Light, dframe2)
summary(ano4)

# Light appears important, possible Light:Regime interaction too, but what about covariates and randeffs...


## glm for covariates

model13 <- glm(SeedCount ~ Light + Regime + Pollinators + fDistance,
              family = poisson
              (link = "log"),
              data = dframe2)

summary(model13)
drop1(model13, test = "Chi")

# Fine, but does pseudoreplication of plot, round and plant have an impact (randeffs)?:


## glmm

model14 <- glmer(SeedCount ~ Light + Regime + Pollinators + Distance
                +(1|fPlot) + (1|fPlantNo) + (1|fRound),
                family = poisson (link = "log"),
                data = dframe2)

summary(model14)
drop1(model14, test = "Chi")

# Obviously first thing to note is that the same problems with confounding variables exist as above, so let's jump straight to the solution:
# first analyse HPS as the major form of lighting

dframe2a <- subset(dframe2,Light!="LED") # HPS + CON in dataset
dframe2a$Regime <- revalue(dframe2a$Regime, c("Control"="None"))
dframe2a$Light  <- revalue(dframe2a$Light, c("CON"="HPS"))

summary(dframe2a)
hist(dframe2a$SeedCount)

# figure within HPS
p16 <- ggplot(dframe2a,aes(x=Regime,y=SeedCount))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d16 <- qplot(Regime, SeedCount, data=dframe2a)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

multiplot(p16, d16, cols=2)

#glmm
model15 <- glmer(SeedCount ~ Regime + Pollinators + Distance
                + (1|fPlantNo) + (1|fRound) + (1|fPlot),
                family = poisson (link = "log"),
                data = dframe2a)

summary(model15)
drop1(model15, test = "Chi") 

# full dataset looked overdispersed, as does histogram, so first check this

mean(dframe2a$SeedCount)
var(dframe2a$SeedCount)

source("OverdispersalFunction.R")
overdisp_fun(model15)

# Yes, overdispersed! Quite badly, so let's try a negative binomial...

model16 <- glmer.nb(SeedCount ~ Regime + Pollinators + Distance
                    + (1|fPlantNo) + (1|fRound) + (1|fPlot),
                    data = dframe2a)

summary(model16)
drop1(model16, test = "Chi")

relgrad <- with(model16@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) # model here does not converge but nothing to worry about

chkres(model16) # good enough...

# No significant effect of Regime on seed count within HPS. What about LED?


dframe2b <- subset(dframe2,Light!="HPS") # LED + CON in dataset
dframe2b$Regime <- revalue(dframe2b$Regime, c("Control"="None"))
dframe2b$Light  <- revalue(dframe2b$Light, c("CON"="LED"))

summary(dframe2b)

# figure within HPS
p17 <- ggplot(dframe2b,aes(x=Regime,y=SeedCount))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d17 <- qplot(Regime, SeedCount, data=dframe2b)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

multiplot(p17, d17, cols=2)

hist(dframe2b$SeedCount)

mean(dframe2b$SeedCount)
var(dframe2b$SeedCount)

model17 <- glmer(SeedCount ~ Regime + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound) + (1|fPlot),
                 family = poisson (link = "log"),
                 data = dframe2b)

summary(model17)
drop1(model17, test = "Chi") 

source("OverdispersalFunction.R")
overdisp_fun(model17)

# Also overdispersed, not quite so bad so let's try QP first...

model18 <- glmmPQL(SeedCount ~   Regime + Pollinators + Distance,
                                     random = list(~1|fPlantNo, ~1|fRound, ~1|fPlot), #Random effects
                                     family = quasipoisson (link = "log"),
                                     data = dframe2b)
summary(model18)
Anova(model18, type = "III")

plot(model18)
chkres(model18) # big patterns in residuals. Probably better try a NB model.


model18a <- glmer.nb(SeedCount ~ Regime + Pollinators + Distance
                     + (1|fPlantNo) + (1|fRound) + (1|fPlot),
                     data = dframe2b)
summary(model18a)
drop1(model18a, test = "Chi")

relgrad <- with(model18a@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) # model here does not converge but nothing to worry about

chkres(model18a) #much better...

# no significant effect of regime on seed count in LED either.
## how about a comparison between LED and HPS?

dframe2c <- subset(dframe2,Regime!="Midnight") # ALL + CON in dataset
dframe2c$Regime <- revalue(dframe2c$Regime, c("Control"="AllNight"))

summary(dframe2c)
hist(dframe2c$SeedCount)

# figure within AllNight
p18 <- ggplot(dframe2c,aes(x=Light,y=SeedCount))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d18 <- qplot(Light, SeedCount, data=dframe2c)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

multiplot(p18, d18, cols=2)

mean(dframe2c$SeedCount)
var(dframe2c$SeedCount) #can tell this is going to be overdispersed!


model19 <- glmer(SeedCount ~ Light + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound) + (1|fPlot),
                 family = poisson (link = "log"),
                 data = dframe2c)

summary(model19)
drop1(model19, test = "Chi") 

overdisp_fun(model19)

# Also overdispersed, not too bad so let's try QP first...

model20 <- glmmPQL(SeedCount ~   Light + Pollinators + Distance,
                   random = list(~1|fPlantNo, ~1|fRound, ~1|fPlot), #Random effects
                   family = quasipoisson (link = "log"),
                   data = dframe2c)
summary(model20)
Anova(model20, type = "III")

plot(model20)
chkres(model20) # again suggestion of a pattern so let's try NB

model20a <- glmer.nb(SeedCount ~ Light + Pollinators + Distance
                     + (1|fPlantNo) + (1|fRound) + (1|fPlot),
                     data = dframe2c)

relgrad <- with(model20a@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) # model here does not converge but nothing to worry about

summary(model20a)
drop1(model20a, test = "Chi")

chkres(model20a) # not amazing but definite improvement

# no effect of light within full-night, finally let's check part-night

dframe2d <- subset(dframe2,Regime!="AllNight") # MID + CON in dataset
dframe2d$Regime <- revalue(dframe2d$Regime, c("Control"="Midnight"))

summary(dframe2d)

# figure within Midnight
p19 <- ggplot(dframe2d,aes(x=Light,y=SeedCount))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d19 <- qplot(Light, SeedCount, data=dframe2d)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

multiplot(p19, d19, cols=2)

hist(dframe2d$SeedCount)
mean(dframe2d$SeedCount)
var(dframe2d$SeedCount) #can tell this is going to be overdispersed too!

model21 <- glmer(SeedCount ~ Light + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound) + (1|fPlot),
                 family = poisson (link = "log"),
                 data = dframe2d)

summary(model21)
drop1(model21, test = "Chi") 

overdisp_fun(model21)

# Also overdispersed, try QP first...

model22 <- glmmPQL(SeedCount ~   Light + Pollinators + Distance,
                   random = list(~1|fPlantNo, ~1|fRound, ~1|fPlot), #Random effects
                   family = quasipoisson (link = "log"),
                   data = dframe2d)
summary(model22)
Anova(model22, type = "III")

# not quite happy with those NaNs, which have occurred because residual DFs are negative (probably too little data?)
# try a NB and see if it's any better...

model23 <- glmer.nb(SeedCount ~ Light + Pollinators + Distance
                    + (1|fPlantNo) + (1|fRound) + (1|fPlot),
                    data = dframe2d)
summary(model23)
drop1(model23, test = "Chi")

relgrad <- with(model23@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) # model here does not converge but nothing to worry about

# chkres function needs debugging for here - in the meantime...
plot(model23)
sresid <- resid(model23, type = "deviance")
hist(sresid)
fitted.glmm <- fitted(model23, level=1)        # Extract the fitted (predicted) values
plot(sresid ~ fitted.glmm)                   # Check for homoscedasticity
binnedplot(fitted(model23),resid(model23))

# these all look fine
# *significant* effect of light within Midnight - seed count significantly higher under HPS than LED

### SO. Only effect of light on seed count is that LED causes less disruption than HPS under midnight switch-offs.


# finally, how about seed weight?

summary(dframe2)

### plots

hist(dframe2$SeedWeight) # check variable distribution - appears to be approx Poisson but not integers so options?
hist(log(dframe2$SeedWeight),10) # Gaussian-ish when logged

dframe2$lSeedWeight <- log(dframe2$SeedWeight,10)
hist(dframe2$lSeedWeight)

# SeedCount

# figure Light
p20 <- ggplot(dframe2,aes(x=Light,y=lSeedWeight))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d20 <- qplot(Light, lSeedWeight, data=dframe2)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

# figure Regime
p21 <- ggplot(dframe2,aes(x=Regime,y=lSeedWeight))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d21 <- qplot(Regime, lSeedWeight, data=dframe2)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

# figure Treatment - Light and Regime combined
p22 <- ggplot(dframe2,aes(x=Treatment,y=lSeedWeight))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d22 <- qplot(Treatment, lSeedWeight, data=dframe2)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

# figure Pollinators
p23 <- ggplot(dframe2,aes(x=Pollinators,y=lSeedWeight))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d23 <- qplot(Pollinators, lSeedWeight, data=dframe2)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

# figure Distance
p24 <- ggplot(dframe2,aes(x=Distance,y=lSeedWeight))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d24 <- qplot(Distance, lSeedWeight, data=dframe2)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

multiplot(p20, d20, p21, d21, p22, d22, p23, d23, p24, d24, cols=5) # no obvious trends, CIs seem to overlap in most cases.


### analysis

## nested anova

ano5 <- aov(lSeedWeight ~ Light*Regime, dframe2)
summary(ano5)

ano6<- aov(lSeedWeight ~ Regime/Light, dframe2)
summary(ano6)

# Nothing really stands out


## glm for covariates

model24 <- glm(lSeedWeight ~ Light + Regime + Pollinators + fDistance,
               family = gaussian
               (link = "identity"),
               data = dframe2)

summary(model24)
drop1(model24, test = "Chi")

# Fine, but does pseudoreplication of plot, round and plant have an impact (randeffs)?:


## glmm

model25 <- lmer(lSeedWeight ~ Light + Regime + Pollinators + Distance
                 +(1|fPlot) + (1|fPlantNo) + (1|fRound),
                 data = dframe2)

summary(model25)
drop1(model25, test = "Chi")

# Obviously first thing to note is that the same problems with confounding variables exist as above, so let's jump straight to the solution:
# first analyse HPS as the major form of lighting

#re-do these to include lSeedWeight
dframe2a <- subset(dframe2,Light!="LED") # HPS + CON in dataset
dframe2a$Regime <- revalue(dframe2a$Regime, c("Control"="None"))
dframe2a$Light  <- revalue(dframe2a$Light, c("CON"="HPS"))

summary(dframe2a)
hist(dframe2a$lSeedWeight)

# figure within HPS
p25 <- ggplot(dframe2a,aes(x=Regime,y=lSeedWeight))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d25 <- qplot(Regime, lSeedWeight, data=dframe2a)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

multiplot(p25, d25, cols=2)

#glmm
model26 <- lmer(lSeedWeight ~ Regime + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound) + (1|fPlot),
                 data = dframe2a)

summary(model26)
drop1(model26, test = "Chi") 

chkres(model26) # OK, though not amazing

# No significant effect of Regime on seed weight within HPS. What about LED?


dframe2b <- subset(dframe2,Light!="HPS") # LED + CON in dataset
dframe2b$Regime <- revalue(dframe2b$Regime, c("Control"="None"))
dframe2b$Light  <- revalue(dframe2b$Light, c("CON"="LED"))

summary(dframe2b)
hist(dframe2b$lSeedWeight)

# figure within HPS
p26 <- ggplot(dframe2b,aes(x=Regime,y=lSeedWeight))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d26 <- qplot(Regime, lSeedWeight, data=dframe2b)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

multiplot(p26, d26, cols=2)

model27 <- lmer(SeedWeight ~ Regime + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound) + (1|fPlot),
                 data = dframe2b)

summary(model27)
drop1(model27, test = "Chi") 

chkres(model27) # Slight trend in sresid vs fitted. Maybe let's have a play with this one...

# remember that lmm(lg(dat)) =/= lg.glmm(dat)

model27a <- glmmPQL(SeedWeight ~ Regime + Pollinators + Distance,
                    random = list(~1|fPlantNo, ~1|fRound, ~1|fPlot),
                    family = gaussian (link = "log"),
                    data = dframe2b)

chkres.PQL(model27a)
#these are better!

summary(model27a)
Anova(model27a, type = "III")

# no significant effect of regime on seed weight in LED either.
## how about a comparison between LED and HPS?

dframe2c <- subset(dframe2,Regime!="Midnight") # ALL + CON in dataset
dframe2c$Regime <- revalue(dframe2c$Regime, c("Control"="AllNight"))

summary(dframe2c)
hist(dframe2c$SeedWeight)
hist(dframe2c$lSeedWeight)

# figure within AllNight
p27 <- ggplot(dframe2c,aes(x=Light,y=lSeedWeight))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d27 <- qplot(Light, lSeedWeight, data=dframe2c)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

multiplot(p27, d27, cols=2)

model28 <- lmer(lSeedWeight ~ Light + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound) + (1|fPlot),
                 data = dframe2c)

summary(model28)
drop1(model28, test = "Chi") 

chkres(model28) # pretty bad. try a PQL

# Also overdispersed, not too bad so let's try QP first...

model28a <- glmmPQL(SeedWeight ~   Light + Pollinators + Distance,
                   random = list(~1|fPlantNo, ~1|fRound, ~1|fPlot), #Random effects
                   family = gaussian (link = "log"),
                   data = dframe2c)

chkres.PQL(model28a) #better

summary(model28a)
Anova(model28a, type = "III")

# no effect of light within full-night, finally let's check part-night

dframe2d <- subset(dframe2,Regime!="AllNight") # MID + CON in dataset
dframe2d$Regime <- revalue(dframe2d$Regime, c("Control"="Midnight"))

summary(dframe2d)

# figure within Midnight
p28 <- ggplot(dframe2d,aes(x=Light,y=lSeedWeight))+
  stat_summary(fun.y="mean",geom="point",alpha=0.7)

d28 <- qplot(Light, lSeedWeight, data=dframe2d)+
  stat_summary(fun.data = "mean_cl_boot", colour = "red")

multiplot(p28, d28, cols=2)

hist(dframe2d$SeedWeight)
hist(dframe2d$lSeedWeight)

model29 <- lmer(lSeedWeight ~ Light + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound) + (1|fPlot),
                 data = dframe2d)

chkres(model29) # again hint of a pattern...

summary(model29)
drop1(model29, test = "Chi") 

model29a <- glmmPQL(SeedWeight ~   Light + Pollinators + Distance,
                   random = list(~1|fPlantNo, ~1|fRound, ~1|fPlot), #Random effects
                   family = gaussian (link = "log"),
                   data = dframe2d)

chkres.PQL(model29a)

summary(model29a)
Anova(model29a, type = "III")

# no effect of light within Midnight

### SO. no effect of light on seed count.

# effect on probability of pollination - driven by regime
# once pollinated, quality is basically conserved (v. limited evidence that HPS is more disruptive than LED from seed counts within midnight regime)
