#Load Libraries
library(corrplot)
library(MASS)
library(effsize)
library(pROC)
library(ISLR)
library(gplots)
library(ROCR)
library(boot)
library(binom)
library(xtable)

#Load in prostate data set
Prostate = read.table(file.choose(), header=T)
dim(Prostate)
sum(is.na(Prostate)) # missing data: 3 in race
Prostate<-na.omit(Prostate); dim(Prostate) # 377 8
Prostate$dcaps = as.factor(Prostate$dcaps)
Prostate$dpros = as.factor(Prostate$dpros)
Prostate=Prostate[ ,-1]
attach(Prostate)
names(Prostate)

#EDA
table(capsule)

#Categorical Data
x1=table(capsule,race)
x2=table(capsule,dpros)
x3=table(capsule,dcaps)
x4=table(capsule,gleason)
x1;x2;x3;x4

#Continuous Data
par(mfrow=c(1,2))
plot(factor(capsule),psa, col="blue", varwidth=T, ylab="PSA Value",xlab="Capsule Penetration")
plot(factor(capsule),age, col="blue", varwidth=T, ylab="Age",xlab="Capsule Penetration")

# d) standardized mean difference: evaluation of relationships between covariates and the response
t.test(age[capsule==0], age[capsule==1])
t.test(psa[capsule==0], psa[capsule==1])
chisq.test(x=race, y=capsule)
chisq.test(x=dpros, y=capsule)
chisq.test(x=dcaps, y=capsule)
chisq.test(x=gleason, y=capsule)

cohen.d(age,factor(capsule))
cohen.d(psa,factor(capsule))


#Interaction Terms
par(mfrow=c(1,3))
plot(factor(gleason),age, col="blue", varwidth=T, ylab="Age",xlab="Gleason Score")
plot(factor(dpros),age, col="blue", varwidth=T, ylab="Age",xlab="Digital Rectal Exam")
plot(factor(race),psa, col="blue", varwidth=T, ylab="PSA",xlab="Race")

par(mfrow=c(1,1))
plot(psa,age)
hist(age)
hist(psa)
hist(log(psa))
plot(log(psa),age)

#Print of Contigency Table
colnames(x1)=cbind("White", "Black")
rownames(x1)=cbind("No", "Yes")
ctrace=xtable(x1,caption="Contigency table comparing race of patient to tumor penetration. 
              Chi-squared analysis determined there is no relationship between tumor penetration and race (p-value = 1).",label="ctrace")
align(ctrace) = "|l|rr|"

colnames(x2)=cbind("None", "Left","Right","Bilobar")
rownames(x2)=cbind("No", "Yes")
ctdpros=xtable(x2,caption="Contigency table compairing the results of digital rectal exam and tumor penetration. 
               Chi-squared analysis determined there is a relationship between exam results and tumor penetration (p-value = 1.974e-08).",label="ctdpros")
align(ctdpros) = "|l|rrrr|"

colnames(x3)=cbind("No", "Yes")
rownames(x3)=cbind("No", "Yes")
ctdcaps=xtable(x3,caption="Contigency table compairing detection of capsular involvement and tumor penetration. 
               Chi-squared analysis determined there is a relationship between capsular involvement and tumor penetration (p-value = 4.221e-06).",label="ctdcaps")
align(ctdcaps) = "|l|rr|"

colnames(x4)=cbind("0","4","5","6","7","8","9")
rownames(x4)=cbind("No", "Yes")
ctgleason=xtable(x4,caption="Contigency table compairing Gleason score and tumor penetration. 
               Chi-squared analysis determined there is a relationship between capsular involvement and tumor penetration (p-value = 2.2e-16).",label="ctgleason")
align(ctgleason) = "|l|rrrrrrr|"

print(ctrace); print(ctdpros); print(ctdcaps); print(ctgleason)


#Model building
#Full Model, then StepAIC to create initial reduced model
FitF<-glm(capsule~.*., family=binomial(link=logit), data=Prostate)

FitR<-stepAIC(FitF)
# b) Parsimonious model: perform backward selection via p-values,

summary(FitF)
summary(FitR)
xtable(print(FitR))

#Parsimonious Model: perform backward selection via p-values

#Candidate model with interaction term
FitRR<-glm(capsule~dpros+psa+gleason+age+age:dpros, family=binomial(link=logit), data=Prostate)
summary(FitRR)

#Candidate model no interaction term
FitRRR<-glm(capsule~dpros+psa+gleason, family=binomial(link=logit),data=Prostate)
summary(FitRRR)

#Compare candidate models to the Full model and each other
anova(FitF,FitRR,test="Chisq")
anova(FitF,FitRRR,test="Chisq")
anova(FitRR,FitRRR,test="Chisq")

#Model evaluation, using examine.logistic.reg function
source("C:/Users/Brian Fury/Stats/Prostate/examine.logistic.reg.R")
mod.fit1<-FitF
examine.logistic.reg(mod.fit1)
mod.fit2<-FitRR
examine.logistic.reg(mod.fit2)
mod.fit3<-FitRRR
examine.logistic.reg(mod.fit3)

# The functions for residual analysis we proposed using are as follows:
one.fourth.root=function(x){
  x^0.25
}

#Interaction Model Outlier Detection
dat.glm <- glm(capsule~dpros+psa+gleason+age+age:dpros, family = binomial, data = Prostate)
dat.mf <- model.frame(dat.glm)
w <- aggregate(formula = capsule~dpros+psa+gleason+age+age:dpros, data = Prostate, FUN = sum)
n <- aggregate(formula = capsule~dpros+psa+gleason+age+age:dpros, data = Prostate, FUN = length)
w.n <- data.frame(w, trials = n$capsule, prop = round(w$capsule/n$capsule,2))
dim(w.n)

# Create EVPs by binning continuous covariates
g = 8 # number of categories
psa_interval = cut(psa, quantile(psa, 0:g/g), include.lowest = TRUE)  # Creates factor with levels 1,2,...,g
levels(psa_interval)

w <- aggregate(formula = capsule ~ psa_interval+gleason+dpros+age+age:dpros, data = Prostate, FUN = sum)
n <- aggregate(formula = capsule ~ psa_interval+gleason+dpros+age+age:dpros, data = Prostate, FUN = length)
w.n <- data.frame(w, trials = n$capsule, prop = round(w$capsule/n$capsule,2))
mod.prelim1 <- glm(formula = capsule/trials ~ psa_interval+gleason+dpros+age+age:dpros,
                   family = binomial(link = logit), data = w.n, weights = trials)
save1=examine.logistic.reg(mod.prelim1, identify.points=T, scale.n=one.fourth.root, scale.cookd=sqrt)
w.n.diag1=data.frame(w.n, pi.hat=round(save1$pi.hat, 2), std.res=round(save1$stand.resid, 2), 
                     cookd=round(save1$cookd, 2), h=round(save1$h, 2))
p=length(mod.prelim1$coef) # number of parameters in model (# coefficients)
ck.out=abs(w.n.diag1$std.res)>2 | w.n.diag1$cookd>4/nrow(w.n) | w.n.diag1$h > 3*p/nrow(w.n)
extract.EVPs=w.n.diag1[ck.out, ]
extract.EVPs

#No interaction model Outlier detection

dat.glm <- glm(capsule~dpros+psa+gleason, family = binomial, data = Prostate)
dat.mf <- model.frame(dat.glm)
w <- aggregate(formula = capsule~dpros+psa+gleason, data = Prostate, FUN = sum)
n <- aggregate(formula = capsule~dpros+psa+gleason, data = Prostate, FUN = length)
w.n <- data.frame(w, trials = n$capsule, prop = round(w$capsule/n$capsule,2))
dim(w.n)

# Create EVPs by binning continuous covariates
g = 8 # number of categories
psa_interval = cut(psa, quantile(psa, 0:g/g), include.lowest = TRUE)  # Creates factor with levels 1,2,...,g
levels(psa_interval)

w <- aggregate(formula = capsule ~ psa_interval+gleason+dpros, data = Prostate, FUN = sum)
n <- aggregate(formula = capsule ~ psa_interval+gleason+dpros, data = Prostate, FUN = length)
w.n <- data.frame(w, trials = n$capsule, prop = round(w$capsule/n$capsule,2))
mod.prelim1 <- glm(formula = capsule/trials ~ psa_interval+gleason+dpros,
                   family = binomial(link = logit), data = w.n, weights = trials)
save1=examine.logistic.reg(mod.prelim1, identify.points=T, scale.n=one.fourth.root, scale.cookd=sqrt)
w.n.diag1=data.frame(w.n, pi.hat=round(save1$pi.hat, 2), std.res=round(save1$stand.resid, 2), 
                     cookd=round(save1$cookd, 2), h=round(save1$h, 2))
p=length(mod.prelim1$coef) # number of parameters in model (# coefficients)
ck.out=abs(w.n.diag1$std.res)>2 | w.n.diag1$cookd>4/nrow(w.n) | w.n.diag1$h > 3*p/nrow(w.n)
extract.EVPs=w.n.diag1[ck.out, ]
extract.EVPs
xtable(extract.EVPs)
#HL Test: Hosmer and Lemeshow goodness-of-fit test using fuction HLtest.R
source("C:/Users/Brian Fury/Stats/Prostate/HLtest.R")
HL = HLTest(FitRR, 10)  # 10 groups by default
# HL test output: Y0 are successes, Y1 are failures
cbind(HL$observed, round(HL$expect, digits=1)) # Observed and Expected table as illustration
HL # HL test results

HL = HLTest(FitRRR, 10)  # 10 groups by default
# HL test output: Y0 are successes, Y1 are failures
cbind(HL$observed, round(HL$expect, digits=1)) # Observed and Expected table as illustration
HL # HL test results

#Interaction Model Odd's Ratio generation

fm2=FitRR 
# construct the table elements
betahat = formatC(signif(fm2$coeff,digits=3), digits=2, format="f", flag="#")
OR = formatC(signif(exp(fm2$coeff),digits=3), digits=2, format="f", flag="#")
SE = formatC(signif(summary(fm2)$coeff[,2],digits=3), digits=2, format="f", flag="#") 
cibounds = formatC(signif(exp(confint(fm2)),digits=3), digits=2, format="f", flag="#") 
pval = formatC(signif(summary(fm2)$coeff[,4],digits=4), digits=4, format="f", flag="#")

x = cbind(betahat, OR, SE, pval, matrix(paste("(", cibounds[,1], ",", cibounds[,2], ")")))
colnames(x) = cbind("Coefficient", "Odds ratio", "SE", "p-value", "95% CI on OR")
rownames(x) = cbind("intercept","dpros", "psa", "gleason", "age*dpros")
inftableFitRR = xtable(x,  
                  caption="Inferences from regressing field goal success on field goal distance,
                  wind conditions, PAT attempt (or not), and if kick resulted in lead change.  The model also 
                  includes interactions of both wind and PAT with distance, respectively. 
                  This inference table was created from the model fit using xtable.", 
                  label="LRinf")
align(inftableFitRR) = "|l|rrrrr|"


#No interaction model Odd's Ratio generation

fm3=FitRRR 
# construct the table elements
betahat = formatC(signif(fm3$coeff,digits=3), digits=2, format="f", flag="#")
OR = formatC(signif(exp(fm3$coeff),digits=3), digits=2, format="f", flag="#")
SE = formatC(signif(summary(fm3)$coeff[,2],digits=3), digits=2, format="f", flag="#") 
cibounds = formatC(signif(exp(confint(fm3)),digits=3), digits=2, format="f", flag="#") 
pval = formatC(signif(summary(fm3)$coeff[,4],digits=4), digits=4, format="f", flag="#")

x = cbind(betahat, OR, SE, pval, matrix(paste("(", cibounds[,1], ",", cibounds[,2], ")")))
colnames(x) = cbind("Coefficient", "Odds ratio", "SE", "p-value", "95% CI on OR")
rownames(x) = cbind("intercept","dpros2","dpros3","dpros4", "psa", "gleason")
inftableFitRRR = xtable(x,  
                  caption="Inferences from regressing field goal success on field goal distance,
                  wind conditions, PAT attempt (or not), and if kick resulted in lead change.  The model also 
                  includes interactions of both wind and PAT with distance, respectively. 
                  This inference table was created from the model fit using xtable.", 
                  label="LRinf")
align(inftableFitRRR) = "|l|rrrrr|"

print(inftableFitRR); print(inftableFitRRR)

#ROC Analysis
#Extract from the fitted model the vector of fitted probabilities
predpr1 <- predict(FitRR,type=c("response"))

predpr2 <- predict(FitRRR,type=c("response"))

predpr3 <- predict(FitF,type=c("response"))

#Generate ROC
roc1<-roc(capsule~predpr1)
roc2<-roc(capsule~predpr2)
roc3<-roc(capsule~predpr3)

#Plot ROC with density Smoothing, generate AUC for actual and smoothed data
par(mfrow=c(1,1))
plot(roc1, main="Interaction Model ROC")
rs1 <- smooth(roc1, method="density")
plot(rs1, add=TRUE, col="blue")
auc(roc1);auc(rs1)

plot(roc2, main="Reduced Model ROC")
rs2 <- smooth(roc2, method="density")
plot(rs2, add=TRUE, col="red")
auc(roc2);auc(rs2)

#Final Comparison of 2 models
roc.test(roc1,roc2)

#Full model ROC, AUC (actual and smoothed), roc test to reduced models
plot(roc3, main="Full Model ROC")
rs3 <- smooth(roc3, method="density")
plot(rs2, add=TRUE, col="purple")
auc(roc3);auc(rs3)

roc.test(roc2,roc3)
roc.test(roc1,roc3)

#Cutoff's, sensitivity, and specificity for all Models; Final Model Delivers highest sensitivity
coords(roc1, "best", ret=c("threshold", "specificity", "sensitivity"))
coords(roc2, "best", ret=c("threshold", "specificity", "sensitivity"))
coords(roc3, "best", ret=c("threshold", "specificity", "sensitivity"))

plot(roc2, print.thres="best", main="Final Model ROC")
rs2 <- smooth(roc2, method="density")
plot(rs2, add=TRUE, col="red")
auc(roc2);auc(rs2)
