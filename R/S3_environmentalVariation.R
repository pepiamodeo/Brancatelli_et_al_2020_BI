
#librerias

library(ggplot2)
library(plyr)
library(multcomp)
library(scales)
library(reshape2)

# load data

DATOS_pp_sitio<-read.csv("./data/data.precipitation.csv",
                    stringsAsFactors = TRUE)
DATOS_pp_sitio$Date<-as.Date(DATOS_pp_sitio$Date)

DATOS_cob_estacion<-read.csv("./data/data.vegetation.csv",
                         stringsAsFactors = TRUE)
DATOS_cob_estacion$Date<-as.Date(DATOS_cob_estacion$Date)

# Environmental variation

## Precipitation

theme_set(theme_bw())
ggplot(data=DATOS_pp_sitio,aes(x=Date,y=PP))+
  geom_point()+
  stat_summary(fun=mean,geom="line")+
  geom_line(aes(group=Site),colour="grey")+
  labs(x= "Date",y="Precipitation (mm)")+
  scale_x_date(labels = date_format("%Y-%m"))

#ggsave("./fig/fig_PP.tiff",width=168,height=168,units="mm",
#       dpi = 600,compression="lzw")

#ggsave("./fig/fig_PP.eps",width=168,height=168,units="mm")

ggsave("./fig/fig_PP.pdf",width=168,height=168,units="mm",
       colormodel = "cmyk")

# Precipitaton values
mean(DATOS_pp_sitio$PP,na.rm = T)
sd(DATOS_pp_sitio$PP,na.rm = T)/length(na.exclude(DATOS_pp_sitio$PP))
range(DATOS_pp_sitio$PP,na.rm = T)

# Termporal correlation between sites
pp_cor<-dcast(DATOS_pp_sitio, Date ~ Site,value.var = "PP")
rownames(pp_cor)<-pp_cor$Date

cor<-cor(na.exclude(pp_cor[,-1]))
cor<-cor[upper.tri(cor)]

mean(cor)
sd(cor)/length(cor)
range(cor)

# Vegetation cover

DATOS_cob_estacion<-transform(DATOS_cob_estacion,
          R_SD=R+SD,
          G_NG=G+NG)
# cover type: A= shrubs, R_SD= rocks and bare soil, G_NG= grasses and herbs

DATOS_cob_estacion<-DATOS_cob_estacion[,-c(5,6,7,8)] # Discard G NG R y SD
# Vegetation cover - orientacion-altitude

DATOS_cob_ggplot<-melt(DATOS_cob_estacion[,-5],id.vars = c("Date","Plot","Site","Site_Plot",
                                                      "Orientation","Altitude"))
ggplot(data=DATOS_cob_ggplot,
       aes(x=Date,y=value,colour=variable))+
    stat_summary(fun = mean,geom="point")+
    stat_summary(fun = mean,geom="line")+
  facet_grid(Altitude~Orientation,scales = "free_y")

theme_set(theme_bw())
ggplot(data=DATOS_cob_ggplot,
       aes(x=Orientation,y=value,fill=variable))+
  stat_summary(fun.y = mean,geom="bar",position="stack")+
  #stat_summary(fun.y = mean_se,geom="errorbar",position="stack")+
  facet_grid(Altitude~.)+
  labs(y="Coverage (%)")+
  scale_fill_grey()

#ggsave("./fig/fig_vegcoverage.tiff",width=174,height=174,units="mm",
#       dpi = 300)
ggsave("./fig/fig_vegcoverage.pdf",width=168,height=168,units="mm",
       colormodel = "cmyk")

# summary of vegetation covers
cob.summary<-ddply(DATOS_cob_estacion,.(Orientation,Altitude),summarise,
          A.mean=mean(A,na.rm=T), # A
          A.SE=sd(A,na.rm=T)/100,
          A.min=min(A,na.rm=T),
          A.max=max(A,na.rm=T),
          G_NG.mean=mean(G_NG,na.rm=T), # G NG
          G_NG.SE=sd(G_NG,na.rm=T)/100,
          G_NG.min=min(G_NG,na.rm=T),
          G_NG.max=max(G_NG,na.rm=T),
          R_SD.mean=mean(R_SD,na.rm=T), # R SD
          R_SD.SE=sd(R_SD,na.rm=T)/100,
          R_SD.min=min(R_SD,na.rm=T),
          R_SD.max=max(R_SD,na.rm=T),
          n=length(R_SD))
cob.summary

# temporal variación in site plot
CV.temporal<-ddply(DATOS_cob_estacion,.(Orientation,Altitude,Site_Plot),summarise,
          A.CV=sd(A,na.rm=T)/mean(A,na.rm=T),
          G_NG.CV=sd(G_NG,na.rm=T)/mean(G_NG,na.rm=T),
          R_SD.CV=sd(R_SD,na.rm=T)/mean(R_SD,na.rm=T))

CV.temporal
# resumen por condición
CV.temporal<-ddply(CV.temporal,.(Orientation,Altitude),summarise,
                   A.CV.mean=mean(A.CV,na.rm=T), # A
                   A.CV.SE=sd(A.CV,na.rm=T)/10,
                   G_NG.CV.mean=mean(G_NG.CV,na.rm=T), # G NG
                   G_NG.CV.SE=sd(G_NG.CV,na.rm=T)/10,
                   R_SD.CV.mean=mean(R_SD.CV,na.rm=T), # R SD
                   R_SD.CV.SE=sd(R_SD.CV,na.rm=T)/10,
                   n=length(R_SD.CV))

summary(c(CV.temporal$A.CV.mean,CV.temporal$G_NG.CV.mean,CV.temporal$R_SD.CV.mean))

#write.csv(CV.temporal,"./tab/CV.cob.temporal.csv")

# Comparison shrubs

DATOS_cob_estacion$Orientation_Altitude<-as.factor(paste(DATOS_cob_estacion$Orientation,DATOS_cob_estacion$Altitude))
DATOS_cob_estacion<-mutate(DATOS_cob_estacion,A_prop = A/100)

fit_cob<-glm(data=DATOS_cob_estacion,A_prop ~ Orientation*Altitude,
             family=binomial)
anova(fit_cob,test="LRT")

fit_cob<-glm(data=DATOS_cob_estacion,A_prop ~ Orientation+Altitude,
             family=binomial)
anova(fit_cob,test="LRT")
summary(fit_cob)


fit_cob<-glm(data=DATOS_cob_estacion,A_prop ~ Orientation,
             family=binomial)
anova(fit_cob,test="LRT")

# Comparison grass and herbs

DATOS_cob_estacion$Orientation_Altitude<-as.factor(paste(DATOS_cob_estacion$Orientation,DATOS_cob_estacion$Altitude))
DATOS_cob_estacion<-mutate(DATOS_cob_estacion,G_NG_prop = G_NG/100)

fit_cob<-glm(data=DATOS_cob_estacion,G_NG_prop ~ Orientation*Altitude,
             family=binomial)
anova(fit_cob,test="LRT")

fit_cob<-glm(data=DATOS_cob_estacion,G_NG_prop ~ Orientation+Altitude,
             family=binomial)
anova(fit_cob,test="LRT")

fit_cob<-glm(data=DATOS_cob_estacion,G_NG_prop ~ Orientation_Altitude,
             family=binomial)
anova(fit_cob,test="LRT")
summary(fit_cob)

summary(glht(fit_cob, mcp(Orientation_Altitude="Tukey")))

comp<-glht(fit_cob, mcp(Orientation_Altitude="Tukey"))
cld(comp)

# Comparison Bare Soil

DATOS_cob_estacion<-mutate(DATOS_cob_estacion,R_SD_prop = R_SD/100)
fit_cob<-glm(data=DATOS_cob_estacion,R_SD_prop ~ Orientation*Altitude,
             family=binomial)
anova(fit_cob,test="LRT")
summary(fit_cob)

fit_cob<-glm(data=DATOS_cob_estacion,R_SD_prop ~ Orientation_Altitude,
             family=binomial)
anova(fit_cob,test="LRT")
summary(fit_cob)

summary(glht(fit_cob, mcp(Orientation_Altitude="Tukey")))
comp<-glht(fit_cob, mcp(Orientation_Altitude="Tukey"))
cld(comp)

# vegetation height

DATOS_cob_estacion_plotHeight<-DATOS_cob_estacion

levels(DATOS_cob_estacion_plotHeight$Orientation)<-c("Northeast Orientation","Southwest Orientation")
levels(DATOS_cob_estacion_plotHeight$Altitude)<-c("Low", "Mid","High")

ggplot(data=DATOS_cob_estacion_plotHeight,
       aes(x=Date,y=veg_height,group=Orientation))+
  stat_summary(fun = mean,geom="point")+
  stat_summary(fun = mean,geom="line",aes(linetype=Orientation))+
  facet_grid(Altitude~.)+
  labs(y="Mean height of vegetation (cm)",x="Date")+
  theme(legend.position = "top")

#ggsave("./fig/fig_vegheight.tiff",width=168,height=168,units="mm",
#       dpi = 600,compression="lzw")
#ggsave("./fig/fig_vegheight.eps", width=168,height=168,units="mm")
ggsave("./fig/fig_vegheight.pdf",width=168,height=168,units="mm",
       colormodel = "cmyk")

ddply(DATOS_cob_estacion,.(Orientation,Altitude),summarise,
                   alt_mean.mean=mean(veg_height,na.rm=T), # A
                   Altm.SE=sd(veg_height,na.rm=T)/100,
                   Altm.min=min(veg_height,na.rm=T),
                   Altm.max=max(veg_height,na.rm=T),
                   n=length(R_SD))

# Comparison Veg Height

fit<-lm(data=DATOS_cob_estacion,veg_height ~ Orientation*Altitude)
anova(fit)
#plot(fit)
summary(fit)

fit<-lm(data=DATOS_cob_estacion,veg_height ~ Orientation_Altitude)
anova(fit)
#plot(fit)
summary(fit)

summary(glht(fit, mcp(Orientation_Altitude="Tukey")))
comp<-glht(fit, mcp(Orientation_Altitude="Tukey"))
cld(comp)

# subset summer
df<-DATOS_cob_estacion[DATOS_cob_estacion$Date>as.Date("2018-12-21")&DATOS_cob_estacion$Date<as.Date("2019-02-15"),]

ggplot(data=df,
       aes(x=Orientation,y=veg_height))+
  stat_summary(fun = mean,geom="bar")+
  facet_grid(Altitude~.)+
  labs(y="Mean height of vegetation (cm)",x="Orientation",
       title="Vegetation height summer 2019")

#ggsave("./fig/fig_vegheight_summer.tiff",width=174,height=174,units="mm",
#       dpi = 300)
ggsave("./fig/fig_vegheight_summer.pdf",width=168,height=168,units="mm",
       colormodel = "cmyk")
