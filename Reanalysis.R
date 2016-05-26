#############################################
#  Binomial logistic regression of seed set #
#############################################

.libPaths(c("C\\rlib", .libPaths()))
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

# also need alternatively coded dataframe

dframeYN <- read.csv("Data\\SeedSetYN.csv")
# dframeYN <- read.csv(file.choose())

names(dframeYN)

dframeYN$fPlot <- factor(dframeYN$Plot) # treat "Plot", "Round", "PlantNo" as factors
dframeYN$fRound <- factor(dframeYN$Round)
dframeYN$fPlantNo <- factor(dframeYN$PlantNo)
dframeYN$fDistance <- factor(dframeYN$Distance) # make "Distance" available as a factor if reqd

dframeYN$Regime <- relevel(dframeYN$Regime,"Control") # relevel variables related to control level
dframeYN$Pollinators <- relevel(dframeYN$Pollinators,"Control")
dframeYN$LitUnlit <- relevel(dframeYN$LitUnlit, "Unlit")

summary(dframeYN)

with(dframeYN,table(Plot,Round)) # seedheads per plot and round - fairly even spread, no obvious patterns, good





library(lme4)  # loading up the libraries
library(car)
library(MASS)
library(praise)
source("CheckResidsFunction.R") # a function for plotting glmer residuals (deviance and Pearson)

### plots

#use install.packages("ggplot2") if necessary
library(ggplot2)
library(scales)
source("MultiplotFunction.R") # function for panel plots in ggplot2 - see http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/


###### Question 1 - complementarity and redundancy of pollination

### Chance of pollination

#  source("MakeRMFunction.R") # a function for reformatting repeated measures - no longer in use in this script

summary(dframe1)

model1 <- lm(cbind(Successes,Failures) ~ Pollinators,
             data = dframe1)

summary(model1)
anova(model1)

# look only within unlit treatment

dframe1c <- subset(dframe1,Light=="CON")
dframe1ca <- subset(dframe1c,Pollinators=="All")
summary(dframe1ca)

dframe1cn <- subset(dframe1c,Pollinators=="Nocturnal")
summary(dframe1cn)

dframe1cd <- subset(dframe1c,Pollinators=="Diurnal")
summary(dframe1cd)

dframe1cc <- subset(dframe1c,Pollinators=="Control")
summary(dframe1cc)


model1c <- lm(cbind(Successes,Failures) ~ Pollinators,
                data = dframe1c)

summary(model1c)
anova(model1c)  # essentially a MANOVA hence the output giving Pillai value

# now, can't do post-hoc tests on repeated measures so let's bring in the YN coded dataset...

summary(dframeYN)

#check it gives similar ANOVA results

modelYN <- lm(SeedSetYN ~ Pollinators,
               data = dframeYN)

summary(modelYN)
anova(modelYN)

library(modelYN)
coefplot(modelYN)

#it does so set up the Tukey test

a1 <- aov(SeedSetYN ~ Pollinators,
           data = dframeYN)
summary(a1)


posthoc <- TukeyHSD(x=a1, 'Pollinators', conf.level=0.95)
posthoc



# repeat for CON only

dframeYNc <- subset(dframeYN,Light=="CON")
summary(dframeYNc)


#check it gives similar ANOVA results

modelYNc <- lm(SeedSetYN ~ Pollinators,
               data = dframeYNc)

summary(modelYNc)
anova(modelYNc)

#it does so set up the Tukey test

a1c <- aov(SeedSetYN ~ Pollinators,
            data = dframeYNc)
summary(a1c)

posthoc <- TukeyHSD(x=a1c, 'Pollinators', conf.level=0.95)
posthoc

praise()

### Seed count & weight

dframe1s <- read.csv("Data\\SeedCount.csv")

names(dframe1s)

dframe1s$fPlot <- factor(dframe1s$Plot) # treat "Plot", "Round", "PlantNo" as factors
dframe1s$fRound <- factor(dframe1s$Round)
dframe1s$fPlantNo <- factor(dframe1s$PlantNo)
dframe1s$fDistance <- factor(dframe1s$Distance) # make "Distance" available as a factor if reqd

dframe1s$Regime <- relevel(dframe1s$Regime,"Control") # relevel variables related to control level
dframe1s$Pollinators <- relevel(dframe1s$Pollinators,"Control")
dframe1s$LitUnlit <- relevel(dframe1s$LitUnlit,"Lit")

summary(dframe1s)

with(dframe1s,table(Plot,Round)) # plants per plot and round - all plants included as they should be

# test

hist(dframe1s$SeedCount) #data approx Poisson, so standard Poisson GLMM may be ok.
plot(dframe1s$SeedCount ~ dframe1s$Pollinators)

model1s <- glmer(SeedCount ~ Pollinators + Distance
                 +(1|fPlantNo) + (1|fRound),
                 family = poisson (link = "log"),
                 data = dframe1s)
summary(model1s)
drop1(model1s, test = "Chi")


#Restrict to unlit
dframe1sc <- subset(dframe1s,Light=="CON")
summary(dframe1sc)

hist(dframe1sc$SeedCount) #data approx Poisson, so standard Poisson GLMM may be ok.
plot(dframe1sc$SeedCount ~ dframe1sc$Pollinators)

model1sc <- glmer(SeedCount ~ Pollinators + Distance
                 +(1|fPlantNo) + (1|fRound),
                 family = poisson (link = "log"),
                 data = dframe1sc)
summary(model1sc)
drop1(model1sc, test = "Chi")

### Seed weight

hist(dframe1s$SeedWeight)
hist(log(dframe1s$SeedWeight),10)

dframe1s$lSeedWeight <- log(dframe1s$SeedWeight,10)
plot(dframe1s$lSeedWeight ~ dframe1s$Pollinators)

model1w <- lmer(lSeedWeight ~ Pollinators + Distance
                +(1|fPlantNo) + (1|fRound),
                data = dframe1s)

summary(model1w)
drop1(model1w, test = "Chi")


# restrict to unlit

hist(dframe1sc$SeedWeight)
hist(log(dframe1sc$SeedWeight),10)

dframe1sc$lSeedWeight <- log(dframe1sc$SeedWeight,10)
plot(dframe1sc$lSeedWeight ~ dframe1sc$Pollinators)

model1wc <- lmer(lSeedWeight ~ Pollinators + Distance
                 +(1|fPlantNo) + (1|fRound),
                 data = dframe1sc)

summary(model1wc)
drop1(model1wc, test = "Chi")


### figures
praise()
#PollinatorsYN

modelYNf <- lmer(SeedSetYN ~ Pollinators + (1|fPlot),
                  data = dframeYN)

summary(modelYNf)
drop1(modelYNf, test = "Chi")



newdata1<-expand.grid(Pollinators=(c("Control","All","Diurnal","Nocturnal")),SeedSetYN=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm1<-model.matrix(terms(modelYNf),newdata1)
newdata1$SeedSetYN = mm1 %*% fixef(modelYNf)
pvar1 <- diag(mm1 %*% tcrossprod(vcov(modelYNf),mm1))
newdata1 <- data.frame(
  newdata1
  , plo = newdata1$SeedSetYN-1.96*sqrt(pvar1)
  , phi = newdata1$SeedSetYN+1.96*sqrt(pvar1)
)

library(plyr)
newdata1$Pollinators <- revalue(newdata1$Pollinators, c("Control"="Caged","All"="Open"))
newdata1  

newdata1$Pollinators<-relevel(newdata1$Pollinators,ref="Open")

#Plot


g1 <- ggplot(newdata1,
              aes(x=Pollinators, y=SeedSetYN, fill=Pollinators))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("gray30","white","gray70","gray50"))+
  scale_y_continuous(labels=percent_format(), limits = c(0,1), oob=squish)+
  guides(fill=FALSE)+
  geom_errorbar(aes(ymin = plo, ymax = phi),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text = element_text(size=19),
        axis.text=element_text(color="black"))+
  labs(x="Pollinator treatment", y="Pollination rate")

g1





###### Question 2 - effect of lighting (lit vs unlit)

model2 <- glmer(cbind(Successes, Failures) ~ LitUnlit + Pollinators + Distance
                + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1)

summary(model2)
drop1(model2, test = "Chi")

chkres(model2)

# nocturnal only

dframe1n <- subset(dframe1,Pollinators=="Nocturnal")

model2n <- glmer(cbind(Successes, Failures) ~ LitUnlit + Distance
                + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1n)

summary(model2n)
drop1(model2n, test = "Chi")

chkres(model2n)

# model failed to converge with fPlot included, but fPlot has very low variance, so removed


# diurnal only

dframe1d <- subset(dframe1,Pollinators=="Diurnal")

model2d <- glmer(cbind(Successes, Failures) ~ LitUnlit + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1d)

summary(model2d)
drop1(model2d, test = "Chi")


# open only

dframe1o <- subset(dframe1,Pollinators=="All")

model2o <- glmer(cbind(Successes, Failures) ~ LitUnlit + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1o)

summary(model2o)
drop1(model2o, test = "Chi")

# control only

dframe1cp <- subset(dframe1,Pollinators=="Control")  #dframe1c already in use

model2c <- glmer(cbind(Successes, Failures) ~ LitUnlit + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1cp)

summary(model2c)
drop1(model2c, test = "Chi")

# Effect of 'being lit' over all variants therein is not significant

### Seed Count

hist(dframe1s$SeedCount)
plot(dframe1s$SeedCount ~ dframe1s$LitUnlit)

model2a <- glmer(SeedCount ~ LitUnlit + Pollinators + Distance
                + (1|fPlantNo) + (1|fRound),
                family = poisson (link = "log"),
                data = dframe1s)

summary(model2a)
drop1(model2a, test = "Chi")

chkres(model2a)

# open pollination

dframe1so <- subset(dframe1s,Pollinators=="All")

hist(dframe1so$SeedCount)
plot(dframe1so$SeedCount ~ dframe1so$LitUnlit)

model2ao <- glmer(SeedCount ~ LitUnlit + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family = poisson (link = "log"),
                 data = dframe1so)

summary(model2ao)
drop1(model2ao, test = "Chi")

chkres(model2ao)


# diurnal pollination

dframe1sd <- subset(dframe1s,Pollinators=="Diurnal")

hist(dframe1sd$SeedCount)
plot(dframe1sd$SeedCount ~ dframe1sd$LitUnlit)

model2ad <- glmer(SeedCount ~ LitUnlit + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1sd)

summary(model2ad)
drop1(model2ad, test = "Chi")

chkres(model2ad)


# nocturnal pollination

dframe1sn <- subset(dframe1s,Pollinators=="Nocturnal")

hist(dframe1sn$SeedCount)
plot(dframe1sn$SeedCount ~ dframe1sn$LitUnlit)

model2an <- glmer(SeedCount ~ LitUnlit + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1sn)

summary(model2an)
drop1(model2an, test = "Chi")

chkres(model2an)


# caged, no pollination

dframe1sc <- subset(dframe1s,Pollinators=="Control")

hist(dframe1sc$SeedCount)
plot(dframe1sc$SeedCount ~ dframe1sc$LitUnlit)

# only 2 data points so no point going further!


### Seed Weight

hist(dframe1s$SeedWeight)
hist(log(dframe1s$SeedWeight,10))

dframe1s$lSeedWeight <- log(dframe1s$SeedWeight,10)

hist(dframe1s$lSeedWeight)
plot(dframe1s$SeedWeight ~ dframe1s$LitUnlit)

model2w <- lmer(lSeedWeight ~ LitUnlit + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound),
                 data = dframe1s)

summary(model2w)
drop1(model2w, test = "Chi")

chkres(model2w)  # dodgy so try a Gamma

model2wa <- glmmPQL(SeedWeight ~ LitUnlit + Pollinators + Distance,
                 random = list(~1|fPlantNo, ~1|fRound),
                 family = quasipoisson (link = "log"),
                 data = dframe1s)

chkres.PQL(model2wa)

summary(model2wa)
Anova(model2wa, type = "III")


# open pollination

dframe1so <- subset(dframe1s,Pollinators=="All")

hist(dframe1so$SeedWeight)
hist(dframe1so$lSeedWeight)
plot(dframe1so$SeedWeight ~ dframe1so$LitUnlit)

model2wo <- glmmPQL(SeedWeight ~ LitUnlit + Distance,
                    random = list(~1|fPlantNo, ~1|fRound),
                    family = quasipoisson (link = "log"),
                    data = dframe1so)

chkres.PQL(model2wo)

summary(model2wo)
Anova(model2wo, type = "III")


# diurnal pollination

dframe1sd <- subset(dframe1s,Pollinators=="Diurnal")

hist(dframe1sd$SeedWeight)
hist(dframe1sd$lSeedWeight)
plot(dframe1sd$SeedWeight ~ dframe1sd$LitUnlit)

model2wd <- glmmPQL(SeedWeight ~ LitUnlit + Distance,
                    random = list(~1|fPlantNo, ~1|fRound),
                    family = quasipoisson (link = "log"),
                    data = dframe1sd)

chkres.PQL(model2wd)

summary(model2wd)
Anova(model2wd, type = "III")


# nocturnal pollination

dframe1sn <- subset(dframe1s,Pollinators=="Nocturnal")

hist(dframe1sn$SeedWeight)
hist(dframe1sn$lSeedWeight)
plot(dframe1sn$SeedWeight ~ dframe1sn$LitUnlit)

model2wn <- glmmPQL(SeedWeight ~ LitUnlit + Distance,
                    random = list(~1|fPlantNo, ~1|fRound),
                    family = quasipoisson (link = "log"),
                    data = dframe1sn)

chkres.PQL(model2wn)

summary(model2wn)
Anova(model2wn, type = "III")


# caged, no pollination

dframe1sc <- subset(dframe1s,Pollinators=="Control")

hist(dframe1sc$SeedWeight)
hist(dframe1sc$lSeedWeight)
plot(dframe1sc$SeedWeight ~ dframe1sc$LitUnlit)

# only 2 data points so no point going further!



### figures

#LitUnlit


#all

model2a <- glmer(SeedCount ~ LitUnlit + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family = poisson (link = "log"),
                 data = dframe1s)

summary(model2a)
drop1(model2a, test = "Chi")

chkres(model2a)

newdata2<-expand.grid(LitUnlit=(c("Lit","Unlit")),Pollinators=(c("Control","All","Diurnal","Nocturnal")),Distance=0,SeedCount=1)
mm2<-model.matrix(terms(model2a),newdata2)
newdata2$SeedCount = mm2 %*% fixef(model2a)
pvar2 <- diag(mm2 %*% tcrossprod(vcov(model2a),mm2))
newdata2 <- data.frame(
  newdata2
  , plo = newdata2$SeedCount-1.96*sqrt(pvar2)
  , phi = newdata2$SeedCount+1.96*sqrt(pvar2)
)
newdata2  

#Plot


g2 <- ggplot(newdata2,
             aes(x=Pollinators, y=SeedCount, fill=LitUnlit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("goldenrod","gray30"))+
  guides(fill=FALSE)+
  xlab("Pollinator treatment")+ ylab("Average seed count per seed capsule")+ 
  geom_errorbar(aes(ymin = plo, ymax = phi),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1))

g2


# pollination treatments separately
#open

model2ao <- glmer(SeedCount ~ LitUnlit + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1so)

summary(model2ao)
drop1(model2ao, test = "Chi")


newdata2ao<-expand.grid(LitUnlit=(c("Lit","Unlit")),Pollinators="All",Distance=0,SeedCount=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2ao<-model.matrix(terms(model2ao),newdata2ao)
newdata2ao$SeedCount = mm2ao %*% fixef(model2ao)
pvar2ao <- diag(mm2ao %*% tcrossprod(vcov(model2ao),mm2ao))
newdata2ao <- data.frame(
  newdata2ao
  , plo = newdata2ao$SeedCount-1.96*sqrt(pvar2ao)
  , phi = newdata2ao$SeedCount+1.96*sqrt(pvar2ao)
)
newdata2ao 

#nocturnal

model2an <- glmer(SeedCount ~ LitUnlit + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1sn)

summary(model2an)
drop1(model2an, test = "Chi")


newdata2an<-expand.grid(LitUnlit=(c("Lit","Unlit")),Pollinators="Nocturnal",Distance=0,SeedCount=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2an<-model.matrix(terms(model2an),newdata2an)
newdata2an$SeedCount = mm2an %*% fixef(model2an)
pvar2an <- diag(mm2an %*% tcrossprod(vcov(model2an),mm2an))
newdata2an <- data.frame(
  newdata2an
  , plo = newdata2an$SeedCount-1.96*sqrt(pvar2an)
  , phi = newdata2an$SeedCount+1.96*sqrt(pvar2an)
)
newdata2an

#diurnal

model2ad <- glmer(SeedCount ~ LitUnlit + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1sd)

summary(model2ad)
drop1(model2ad, test = "Chi")


newdata2ad<-expand.grid(LitUnlit=(c("Lit","Unlit")),Pollinators="Diurnal",Distance=0,SeedCount=1) #Sum.Of.Count doesn't affect predict function and is needed for model.matrix
mm2ad<-model.matrix(terms(model2ad),newdata2ad)
newdata2ad$SeedCount = mm2ad %*% fixef(model2ad)
pvar2ad <- diag(mm2ad %*% tcrossprod(vcov(model2ad),mm2ad))
newdata2ad <- data.frame(
  newdata2ad
  , plo = newdata2ad$SeedCount-1.96*sqrt(pvar2ad)
  , phi = newdata2ad$SeedCount+1.96*sqrt(pvar2ad)
)
newdata2ad


# stitch together

newdata2a <- rbind(newdata2ao,newdata2ad,newdata2an)
newdata2a

newdata2a$Pollinators <- revalue(newdata2a$Pollinators, c("All"="Open"))
newdata2a  


#Plot


g2a <- ggplot(newdata2a,
             aes(x=Pollinators, y=SeedCount, fill=LitUnlit))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("white","gray30"))+
  scale_y_continuous(limits = c(0,6), breaks = c(0,1,2,3,4,5,6))+
  xlab("Pollinator treatment")+ ylab("Average log seed count per seed capsule")+ 
  geom_errorbar(aes(ymin = plo, ymax = phi),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g2a






###### Question 3 - effect of lighting treatments (LED vs HPS & FN vs PN)

### chance of pollination

# overall

dframe1a <- subset(dframe1,Regime!="Control") # Midnight + AllNight
summary(dframe1a)

model3 <- glmer(cbind(Successes, Failures) ~ Light + Regime + Pollinators + Distance
                + (1|fPlantNo) + (1|fRound),
                family = binomial (link = "logit"),
                data = dframe1a)

summary(model3)
drop1(model3, test = "Chi")

chkres(model3) # fine


# open only

dframe1ao <- subset(dframe1a,Pollinators=="All")

model3o <- glmer(cbind(Successes, Failures) ~ Light + Regime + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1ao)

summary(model3o)
drop1(model3o, test = "Chi")



# diurnal only

dframe1ad <- subset(dframe1a,Pollinators=="Diurnal")

model3d <- glmer(cbind(Successes, Failures) ~ Light + Regime + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1ad)

summary(model3d)
drop1(model3d, test = "Chi")


# nocturnal only

dframe1an <- subset(dframe1a,Pollinators=="Nocturnal")

model3n <- glmer(cbind(Successes, Failures) ~ Light + Regime + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1an)

summary(model3n)
drop1(model3n, test = "Chi")

chkres(model3n)


# control only

dframe1ac <- subset(dframe1a,Pollinators=="Control")  

model3c <- glmer(cbind(Successes, Failures) ~ Light + Regime + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family = binomial (link = "logit"),
                 data = dframe1ac)

summary(model3c)
drop1(model3c, test = "Chi")


### Seed Count

dframe1as <- subset(dframe1s,Regime!="Control") # Midnight + AllNight
summary(dframe1as)

hist(dframe1as$SeedCount)
plot(dframe1as$SeedCount ~ dframe1as$Light)
plot(dframe1as$SeedCount ~ dframe1as$Regime)

model3a <- glmer(SeedCount ~ Light + Regime + Pollinators + Distance
                 + (1|fPlantNo) + (1|fRound),
                 family = poisson (link = "log"),
                 data = dframe1as)

summary(model3a)
drop1(model3a, test = "Chi")

chkres(model3a)  # possible overdispersion - try quasi-poisson

source("OverdispersalFunction.R")
overdisp_fun(model3a)


# open pollination

dframe1aso <- subset(dframe1as,Pollinators=="All")

hist(dframe1aso$SeedCount)
plot(dframe1aso$SeedCount ~ dframe1aso$Light)
plot(dframe1aso$SeedCount ~ dframe1aso$Regime)

model3ao <- glmer(SeedCount ~ Light + Regime + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1aso)

summary(model3ao)
drop1(model3ao, test = "Chi")

chkres(model3ao)


# diurnal pollination

dframe1asd <- subset(dframe1as,Pollinators=="Diurnal")

hist(dframe1asd$SeedCount)
plot(dframe1asd$SeedCount ~ dframe1asd$Light)
plot(dframe1asd$SeedCount ~ dframe1asd$Regime)

model3ad <- glmer(SeedCount ~ Light + Regime + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1asd)

summary(model3ad)
drop1(model3ad, test = "Chi")

chkres(model3ad)

# nocturnal pollination

dframe1asn <- subset(dframe1as,Pollinators=="Nocturnal")

hist(dframe1asn$SeedCount)
plot(dframe1asn$SeedCount ~ dframe1asn$Light)
plot(dframe1asn$SeedCount ~ dframe1asn$Regime)

model3an <- glmer(SeedCount ~ Light + Regime + Distance
                  + (1|fPlantNo) + (1|fRound),
                  family = poisson (link = "log"),
                  data = dframe1asn)

summary(model3an)
drop1(model3an, test = "Chi")

chkres(model3an)


# caged, no pollination

dframe1asc <- subset(dframe1as,Pollinators=="Control")

hist(dframe1asc$SeedCount)
plot(dframe1asc$SeedCount ~ dframe1asc$Light)
plot(dframe1asc$SeedCount ~ dframe1asc$Regime)

# only 2 data points so no point going further!


### Seed Weight

hist(dframe1as$SeedWeight)
hist(log(dframe1as$SeedWeight,10))

dframe1as$lSeedWeight <- log(dframe1as$SeedWeight,10)

hist(dframe1as$lSeedWeight)
plot(dframe1as$SeedWeight ~ dframe1as$Light)
plot(dframe1as$SeedWeight ~ dframe1as$Regime)

model3w <- lmer(lSeedWeight ~ Light + Regime + Pollinators + Distance
                + (1|fPlantNo) + (1|fRound),
                data = dframe1as)

summary(model3w)
drop1(model3w, test = "Chi")

chkres(model3w)  


# open pollination

dframe1aso <- subset(dframe1as,Pollinators=="All")

hist(dframe1aso$SeedWeight)
hist(dframe1aso$lSeedWeight)
plot(dframe1aso$SeedWeight ~ dframe1aso$Light)
plot(dframe1aso$SeedWeight ~ dframe1aso$Regime)

model3wo <- lmer(lSeedWeight ~ Light + Regime + Distance
                 + (1|fPlantNo) + (1|fRound),
                 data = dframe1aso)

chkres(model3wo)

summary(model3wo)
drop1(model3wo, test = "Chi")


# diurnal pollination

dframe1asd <- subset(dframe1as,Pollinators=="Diurnal")

hist(dframe1asd$SeedWeight)
hist(dframe1asd$lSeedWeight)
plot(dframe1asd$SeedWeight ~ dframe1asd$Light)
plot(dframe1asd$SeedWeight ~ dframe1asd$Regime)

model3wd <- lmer(lSeedWeight ~ Light + Regime + Distance
                 + (1|fPlantNo) + (1|fRound),
                 data = dframe1asd)

chkres(model3wd)

summary(model3wd)
drop1(model3wd, test = "Chi")


# nocturnal pollination

dframe1asn <- subset(dframe1as,Pollinators=="Nocturnal")

hist(dframe1asn$SeedWeight)
hist(dframe1asn$lSeedWeight)
plot(dframe1asn$SeedWeight ~ dframe1asn$Light)
plot(dframe1asn$SeedWeight ~ dframe1asn$Regime)

model3wn <- lmer(lSeedWeight ~ Light + Regime + Distance
                 + (1|fPlantNo) + (1|fRound),
                 data = dframe1asn)

chkres(model3wn)

summary(model3wn)
drop1(model3wn, test = "Chi")

# caged, no pollination

dframe1asc <- subset(dframe1as,Pollinators=="Control")

hist(dframe1asc$SeedWeight)
hist(dframe1asc$lSeedWeight)
plot(dframe1asc$SeedWeight ~ dframe1asc$Light)
plot(dframe1asc$SeedWeight ~ dframe1asc$Regime)

# only 2 data points so no point going further!




### figures

#Light + Regime

#all

# need alternative coding of data

dframeYNa <- subset(dframeYN,Regime!="Control") # Midnight + AllNight
summary(dframeYNa)

dframeYNa$Regime <- relevel(dframeYNa$Regime,"Midnight")
dframeYNa$Light <- relevel(dframeYNa$Light,"LED")


#revised model

model3YN <- glmer(SeedSetYN ~ Light + Regime + Distance
                   + (1|fPlantNo) + (1|fRound),
                   family = binomial (link = "logit"),
                   data = dframeYNa)

summary(model3YN)
drop1(model3YN, test = "Chi")   # very similar outputs so fine to proceed



newdata3<-expand.grid(Light=(c("HPS","LED")),Regime=(c("AllNight","Midnight")),Distance=0,SeedSetYN=1)
mm3<-model.matrix(terms(model3YN),newdata3)
newdata3$SeedSetYN = exp(mm3 %*% fixef(model3YN))
pvar3 <- diag(mm3 %*% tcrossprod(vcov(model3YN),mm3))
newdata3 <- data.frame(
  newdata3
  , plo = newdata3$SeedSetYN-1.96*sqrt(pvar3)
  , phi = newdata3$SeedSetYN+1.96*sqrt(pvar3)
)
newdata3 



#Plot


g3 <- ggplot(newdata3,
              aes(x=Light, y=SeedSetYN, fill=Regime))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("goldenrod","gray30"))+
  guides(fill=FALSE)+
  xlab("Pollinator treatment")+ ylab("Average seed count per seed capsule")+ 
  geom_errorbar(aes(ymin = plo, ymax = phi),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1))

g3


#alternate method

model3YNs <- glm(SeedSetYN ~ Light + Regime + Distance,
                 family = binomial (link = logit),
                 data = dframeYNa)

summary(model3YNs)
drop1(model3YNs, test = "Chi")   # very similar outputs so fine to proceed

newdata3a<-expand.grid(Light=(c("HPS","LED")),Regime=(c("AllNight","Midnight")),Distance=0,SeedSetYN=1)
newdata3a$SeedSetYN <- predict(model3YNs,newdata=newdata3a,type="response")
preddat <- predict(model3YNs,newdata=newdata3a,type="response",se.fit=TRUE)
preddat

newdata3a <- cbind(newdata3a,preddat)
newdata3a

newdata3a <- data.frame(
  newdata3a
  , plo = newdata3a$fit-1.96*newdata3a$se.fit
  , phi = newdata3a$fit+1.96*newdata3a$se.fit
)
newdata3a

newdata3a$Regime <- revalue(newdata3a$Regime, c("AllNight"="Full night","Midnight"="Part night"))
newdata3a

#Plot


g3a <- ggplot(newdata3a,
             aes(x=Regime, y=SeedSetYN, fill=Light))+
  geom_bar(colour="black",stat="identity",position=position_dodge())+
  scale_fill_manual(values=c("gray70","gray30"))+
  scale_y_continuous(labels=percent_format(), limits = c(0,1), oob=squish)+
  xlab("Lighting regime")+ ylab("Pollination rate")+ 
  geom_errorbar(aes(ymin = plo, ymax = phi),
                colour="black",size=0.8, width=0.5,
                position=position_dodge(width=0.9))+
  theme(panel.background=element_rect(fill="white"),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(colour="gray70"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color="black",fill=F,size=1),
        text=element_text(size=19),
        axis.text=element_text(color="black"),
        legend.title=element_blank())

g3a

