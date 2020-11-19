
# Time to event analysis - Seedling Emergence

# libraries
library(survival)
library(ggplot2)
library(survminer)
library(HH)

# load data

DATOS.EME<-read.csv("./data/data.emergence.csv",
                    stringsAsFactors = TRUE)
DATOS.EME$Date<-as.Date(DATOS.EME$Date)

# ASSUMPTIONS

# Define function mlogmlog() to calculate -log(-log(S(t)))
mlogmlog <- function(y){-log(-log(y))}

# Estimate Kaplan-Meier survivor functions for each of the levels of one factor

fit.campo <- survfit(Surv(Days,Emergence) ~ Orientation+Altitude, type="kaplan-meier",
                    data=DATOS.EME)
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

# Fit a generalized linear model predicting days from site and date
multicollinearitycheck <- glm(Days ~ Orientation*Altitude, data=DATOS.EME)
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

# Load survival library

# SINGLE MODELS

## ORIENTATION
cox.campo <- coxph(Surv(Days, Emergence) ~ Orientation+frailty(Site_Plot), data=DATOS.EME)
print(summary(cox.campo))

km.fit<-survfit(Surv(Days, Emergence) ~ Orientation, 
                DATOS.EME, type="kaplan-meier") 

# Plot the survivor function
y_inv<-function(y){1-y}
plot(km.fit, fun=y_inv, conf.int=F, xlab="time [days]",
     ylab="probability of germinating",mark.time = F,
     col=c("red","blue"))

## ALTITUDE

cox.campo <- coxph(Surv(Days, Emergence) ~ Altitude+frailty(Site_Plot), data=DATOS.EME)
print(summary(cox.campo))


km.fit<-survfit(Surv(Days, Emergence) ~ Altitude, 
                DATOS.EME, type="kaplan-meier") 

# Plot the survivor function
y_inv<-function(y){1-y}
plot(km.fit, fun=y_inv, conf.int=F, xlab="time [days]",
     ylab="probability of germinating",mark.time = F,
     lty=c("solid","longdash","dotted"))

# MODEL WITH TWO EXPLANATORY VARIABLES

## ORIENTATION AND ALTITUDE

cox.campo <- coxph(Surv(Days, Emergence) ~ Orientation*Altitude+frailty(Site_Plot), data=DATOS.EME)
print(summary(cox.campo))


km.fit<-survfit(Surv(Days, Emergence) ~ Orientation+Altitude, 
                DATOS.EME, type="kaplan-meier") 

# Plot the survivor function
y_inv<-function(y){1-y}

p<-ggsurvplot(
        km.fit,  
        data = DATOS.EME,
        fun = y_inv,
        conf.int = TRUE,         # show confidence intervals for point estimates of survival curves.
        xlab = "Time (days)",
        ylab = "Cumulative Emergence Probability",
        ggtheme = theme_bw(),
        break.time.by = 100,
        legend.labs = c("Low", "Mid","High","Low", "Mid","High"))+
        guides(fill = guide_legend(title = 'Altitude'),
               color = guide_legend(title = 'Altitude'))

# Labels for Facet Orientation
labs<-c("Northeast Orientation","Southwest Orientation")
names(labs)<-c("N","S")

p_eme<-p$plot+facet_wrap(~Orientation,labeller = labeller(Orientation=labs))

theme_set(theme_bw())
print(p_eme)

# output plots in different formats

#ggsave("./fig/fig_emergence.tiff",width=168,height=130,units="mm",
#       dpi = 600,compression="lzw")

#ggsave("./fig/fig_emergence.eps",width=168,height=130,units="mm")

ggsave("./fig/fig_emergence.pdf",width=168,height=130,units="mm",
       colormodel = "cmyk")
