
# Time to event analysis - Seedling survival

# library
library(survival)
library(ggplot2)
library(survminer)
library(HH)

# load data

DATOS.SUP<-read.csv("./data/data.survival.csv",
                    stringsAsFactors = TRUE)
DATOS.SUP$Date<-as.Date(DATOS.SUP$Date)

# SUPUESTOS

# Define function mlogmlog() to calculate -log(-log(S(t)))
mlogmlog <- function(y){-log(-log(y))}

# Estimate Kaplan-Meier survivor functions for each of the levels of one factor

fit.campo <- survfit(Surv(Days.mort,Mortality) ~ Orientation+Altitude, type="kaplan-meier",
                     data=DATOS.SUP)
summary(fit.campo)
# Plot -log(-log(S(t))) versus log(t)

plot(fit.campo, mark.time=F,fun=mlogmlog, log="x",
     xlab="t [days]", ylab="-log(-log(S(t)))",
     lty=c("solid","longdash","dotted"), 
     col=c("red","red","red","blue","blue","blue"),lwd=1.75)

# McNair et al 2013. Different covariate values will produce functions with the same shape but 
# different elevations. The Cox model is reasonably robust, 
# so it is only necessary to worry about clear departures from proportionality, 
# as indicated by decisive crossing of the functions for 
# two or more covariates in the diagnostic plot

###### Multicollinearity ####################################
# variance inflation factor, similar to standard multiple regression analysis

# Load HH library

# Read in data

# Fit a generalized linear model predicting days from site and date
multicollinearitycheck <- glm(Days.mort ~ Orientation*Altitude, data=DATOS.SUP)
# Check for multicollinearity among covariates
vif <- vif(multicollinearitycheck)
# Print the results on the screen
print(vif)

# Rhelp:vif. 
# The VIF for predictor i is 1/(1-R_i^2), where R_i^2 is the R^2 from a 
# regression of predictor i against the remaining predictors. 
# If R_i^2 is close to 1, this means that predictor i is well explained by a 
# linear function of the remaining predictors, and, therefore, 
# the presence of predictor i in the model is redundant. 
# Values of VIF exceeding 5 are considered evidence of collinearity: 
# The information carried by a predictor having such a VIF is contained 
# in a subset of the remaining predictors. 
# If, however, all of a model's regression coefficients differ significantly 
# from 0 (p-value < .05), a somewhat larger VIF may be tolerable.

############## Cox model building ###############

# SINGLE MODELS

## ORIENTATION
cox.campo <- coxph(Surv(Days.mort, Mortality) ~ Orientation+frailty(Site_Plot), data=DATOS.SUP)
print(summary(cox.campo))

km.fit<-survfit(Surv(Days.mort, Mortality) ~ Orientation, 
                DATOS.SUP, type="kaplan-meier") 

# Plot the survivor function

plot(km.fit, conf.int=F, xlab="time [days]",
     ylab="probability of survival",mark.time = F,
     col=c("red","blue"))

## ALTITUDE

cox.campo <- coxph(Surv(Days.mort, Mortality) ~ Altitude+frailty(Site_Plot), data=DATOS.SUP)
print(summary(cox.campo))


km.fit<-survfit(Surv(Days.mort, Mortality) ~ Altitude, 
                DATOS.SUP, type="kaplan-meier") 

# Plot the survivor function

plot(km.fit, conf.int=F, xlab="time [days]",
     ylab="probability of survival",mark.time = F,
     lty=c("solid","longdash","dotted"))

# MODEL WITH TWO VARIABLES

## ORIENTATION AND ALTITUDE

cox.campo <- coxph(Surv(Days.mort, Mortality) ~ Orientation*Altitude+frailty(Site_Plot), data=DATOS.SUP)
print(summary(cox.campo))

# Number of survived
cox.campo$n-cox.campo$nevent
# Overall survived proportion
1-(cox.campo$nevent/cox.campo$n)
# Days from sown (total experiment)
summary(DATOS.SUP$Days)
# Days from emergence
summary(DATOS.SUP$Days.mort)

km.fit<-survfit(Surv(Days.mort, Mortality) ~ Orientation+Altitude, 
                DATOS.SUP, type="kaplan-meier") 

km.fit
summary(km.fit)

# Survival at one year

fit365<-survfit(Surv(Days.mort, Mortality) ~ 1, 
                DATOS.SUP, type="kaplan-meier") 

fit365
summary(fit365)

# Dates of median days of mortality
DATOS.SUP[DATOS.SUP$Days.mort%in%c(49,190,200,106,64,141,106),"Date"]

days.mort<-rep(x=km.fit$time,times=km.fit$n.event)
summary(days.mort)

dates.mort<-DATOS.SUP[DATOS.SUP$Days.mort%in%days.mort,"Date"]
summary(dates.mort)

sum(round(table(dates.mort)/cox.campo$nevent,digits = 2)[22:27])

grid.arrange(
ggplot(as.data.frame(dates.mort),aes(x=dates.mort))+
        geom_histogram(bins = 18)+
        labs(x="Date"),
ggplot(as.data.frame(days.mort),aes(x=days.mort))+
        geom_histogram()+
        labs(x="Days since emergence")
)

# Plot the survivor function

p<-ggsurvplot(
        km.fit,   
        data = DATOS.SUP,
        conf.int = TRUE,         # show confidence intervals for point estimates of survival curves.
        xlab = "Time (days)",
        ggtheme = theme_bw(),
        surv.median.line = "v",  # add the median survival pointer.
        break.time.by = 100,
        legend.labs = c("Low", "Mid","High","Low", "Mid","High"))+
        guides(fill = guide_legend(title = 'Altitude'),
               color = guide_legend(title = 'Altitude'))

# Labels for Facet Orientation
labs<-c("Northeast Orientation","Southwest Orientation")
names(labs)<-c("N","S")

p_sup<-p$plot+facet_wrap(~Orientation,labeller = labeller(Orientation=labs))

theme_set(theme_bw())
print(p_sup)

#ggsave("./fig/fig_survival.tiff",width=168,height=130,units="mm",
#       dpi = 600,compression="lzw")

#ggsave("./fig/fig_survival.eps",width=168,height=130,units="mm")

ggsave("./fig/fig_survival.pdf",width=168,height=130,units="mm",
       colormodel = "cmyk")
