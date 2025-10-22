rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

libs<-c("ggplot2","dplyr","ggfortify","fitdistrplus","brms","future","car", "GGally","gamlss.dist",
            "ggeffects","rstan","ggbreak", "cowplot", "doBy", "DHARMa","qqplotr", "bayesplot","arm", "psych", "posterior")

lapply(libs,FUN=require,character.only=T)

# Analisis de datos de tiempo hasta el evento
# Los datos corresponden a tiempo hasta la mortalidad en moscas a las que se es aplicó insecticida
# Se aplico sobre 4 lineas de moscas (LINE 1 a 4) halladas como resistentes o susceptibles en un experimento anterior
# El factor tratamiento es la aplicacion o no de un posible sinergista al insecticida


data<-read.table("data.csv", sep=",", header=T)
str(data)
data$LINE<-as.factor(data$LINE)


formw=bf(TIME | cens(censor)  ~ Treatment*SEX*LINE+(1|ID)+s(TIME,k=5, by=Treatment),family="weibull")

get_prior(formw,data=data , family="weibull")

prior.m2<-c(set_prior("gamma(log(10), 1)",class = "shape"),#Default priors for shape are too far to the structure of this data, causes convergenge problems
            set_prior("normal(0, 2)",class = "b"))

########## ANALISIS ############

set.seed(123)
rstan::rstan_options(disable_march_warning = TRUE)#
data$LINE<-as.factor(data$LINE)

options(mc.cores = parallel::detectCores())

str(data)

mw<-datamw<- brm(data = data, 
          family='weibull',
          prior=prior.m2,
          TIME | cens(censor)  ~ Treatment*LINE*SEX+(1|ID),
          warmup = 1000,chains=3,thin=3, # Aumente thin porque habia autocorrelacion en los parametros de efectos fijos
         iter=6000, refresh = 0,control=list(adapt_delta= 0.95))# Tuve que aumentar adapt_delta y max_treedepth porque habia divergent transitions

######## VALIDACION ####
summary(mw)

neff_ratio(mw)[neff_ratio(mw)<0.7]

windows()
plot(mw, variable = "shape")

nuts <- nuts_params(mw)
draws_array <- as_draws_array(mw)

windows()
mcmc_scatter(
  as_draws_array(mw),
  pars = c("b_Intercept", "b_TreatmentP"),
  np = nuts
)


posterior <- as_draws_df(mw,ndraws=10)

windows()
pairs(posterior[, c("b_Intercept", "b_TreatmentP", "b_LINE3", "b_LINE2", "b_TreatmentP:LINE3", "b_TreatmentP:LINE2", "shape", "lprior",  "lp__")])


#Convergencia de las cadenas
v<-variables(mw)[grep("r_",variables(mw))]

windows()
mcmc_dens_overlay(mw, regex_pars = c("b"))+
  geom_density(lwd=1.2, alpha=0.1)+ 
  theme_bw()+
  ylab("Probability density")+
  theme(legend.position="none",
        axis.text=element_text(size = 14),  
        axis.title = element_text(size = 18),
        strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size = 12)) 


windows()
mcmc_dens_overlay(mw, pars = v[1:36])+#HAY que usar pars y no regex_pars sino no las filtra
  geom_density(lwd=1.2, alpha=0.1)+ 
  theme_bw()+
  ylab("Probability density")+
  theme(legend.position="none",
        axis.text=element_text(size = 14),  
        axis.title = element_text(size = 18),
        strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size = 12)) 

windows()
mcmc_dens_overlay(mw, pars = v[100:120])+#HAY que usar pars y no regex_pars sino no las filtra
  geom_density(lwd=1.2, alpha=0.1)+ 
  theme_bw()+
  ylab("Probability density")+
  theme(legend.position="none",
        axis.text=element_text(size = 14),  
        axis.title = element_text(size = 18),
        strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size = 12)) 


windows()
mcmc_dens_overlay(mw, pars = v[200:220])+#HAY que usar pars y no regex_pars sino no las filtra
  geom_density(lwd=1.2, alpha=0.1)+ 
  theme_bw()+
  ylab("Probability density")+
  theme(legend.position="none",
        axis.text=element_text(size = 14),  
        axis.title = element_text(size = 18),
        strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size = 12)) 

# Autocorrelacion de los valores muestrados de los parametros por cadena
#Si hay autocorrelacion, es necesario aumentar el thining

windows()
mcmc_acf(mw, pars = v[1:12])+
  geom_line(size=0.8)+
  scale_x_continuous(limits=c(1,10), breaks=seq(1,9,2))+
  scale_y_continuous(limits=c(-0.2,0.6))+
  theme(axis.text.y=element_text(size = 16),
        axis.text.x=element_text(size = 14, angle=90),
        axis.title = element_text(size = 18))

windows()
mcmc_acf(mw, pars = v[100:110])+
  geom_line(size=0.8)+
  scale_x_continuous(limits=c(1,10), breaks=seq(1,9,2))+
  scale_y_continuous(limits=c(-0.2,0.6))+
  theme(axis.text.y=element_text(size = 16),
        axis.text.x=element_text(size = 14, angle=90),
        axis.title = element_text(size = 18))

windows()
mcmc_acf(mw, regex_pars = c("b"))+
  geom_line(size=0.8)+
  scale_x_continuous(limits=c(1,10), breaks=seq(1,9,2))+
  scale_y_continuous(limits=c(-0.2,0.6))+
  theme(axis.text.y=element_text(size = 16),
        axis.text.x=element_text(size = 14, angle=90),
        axis.title = element_text(size = 18))


## Mixing de las cadenas

windows()
mcmc_trace(mw, pars = v[1:12])+
  theme_bw()+
  ylab("Probability density")+
  theme(legend.position="none",
        axis.text=element_text(size = 14),  
        axis.title = element_text(size = 18),
        strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size = 12)) 

windows()
mcmc_trace(mw, regex_pars = c("b"))+
  theme_bw()+
  ylab("Probability density")+
  theme(legend.position="none",
        axis.text=element_text(size = 14),  
        axis.title = element_text(size = 18),
        strip.background=element_rect(colour = "black", fill = "white"),
        strip.text = element_text(size = 12)) 

########### PP-checks ###########

mw$data$event<-ifelse(mw$data$censor=="right", "0", "1")
mw$data$event<-as.numeric(mw$data$event)

mw$data$TL<-paste(mw$data$Treatment,mw$data$LINE)

yrep <- posterior_predict(mw,ndraws=50)

windows()
ppc_km_overlay_grouped(y=mw$data$TIME,yrep,status_y=mw$data$event, group = mw$data$LINE)

windows()
ppc_km_overlay_grouped(y=mw$data$TIME,yrep,status_y=mw$data$event, group = mw$data$Treatment)

windows()
ppc_km_overlay_grouped(y=mw$data$TIME,yrep,status_y=mw$data$event, group = mw$data$SEX)

windows()
ppc_km_overlay_grouped(y=mw$data$TIME,yrep,status_y=mw$data$event, group = mw$data$TL)

###### Visualización de resultados #######

a<-conditional_effects(mw, prob=0.9 ,conditions=make_conditions(mw$data, "SEX", "Treatment", "LINE") )

a<-rbind(a$`Treatment:LINE`,a$`Treatment:SEX`,a$`LINE:SEX`)

a$LT<-paste(a$LINE,a$Treatment)

a$LTS<-paste(a$LINE,a$Treatment,a$SEX)

a <- a %>% distinct(LTS, .keep_all = TRUE)


a$Treatment<-as.character(a$Treatment)
a$Treatment[a$Treatment=="P"]<-"PBO"
a$Treatment[a$Treatment=="C"]<-"Control"

levels(a$SEX)<-c("Hembras", "Machos")

a$RES<-"Resistente"

a$RES[a$LINE=="1" | a$LINE=="100"]<-"Susceptible"


windows()
ggplot(a, aes(x=LINE, y=estimate__, color=SEX,group=LTS, shape=Treatment)) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__),size = 2,position = "dodge2", 
                width=.2) +
  geom_point(position = position_dodge2(width=.2),size = 5)+
  ylab("TIME (h)") +
  xlab("") +
  ggtitle("")+
  guides(color=guide_legend(title="SEX"),shape=guide_legend(title="Treatmentamiento"))+
  theme_classic(base_size = 25)+
  theme(legend.position="bottom")


windows()
ggplot(a, aes(x=Treatment, y=estimate__, color=LINE, group=LTS, shape=RES)) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__),size = 3,position = "dodge2", 
                width=.2) +
  geom_point(position = position_dodge2(width=.2),size = 6)+
  facet_grid(cols=vars(SEX)) +
  ylab("TIME (h)") +
  xlab("") +
  ggtitle("Spinosad")+
  guides(color=guide_legend(title="Línea"), shape=guide_legend(title=""))+
  theme_minimal(base_size = 20)+
  theme(legend.position = "bottom")

windows()
conditional_effects(mw, effects = "LINE", prob= 0.95 )

windows()
conditional_effects(mw, effects = "Treatment" )

windows()
conditional_effects(mw, effects = "SEX")

####### CURVAS #########

newdata <- expand.grid(
  Treatment = unique(mw$data$Treatment),
  LINE = unique(mw$data$LINE),
  SEX = unique(mw$data$SEX)
)

# Predictor LINEl
n_draws <- 100

# 2) extraer matriz [draws × nrow(newdata)] de λ
lambda_draws <- posterior_linpred(
  mw,
  newdata    = newdata,
  transform  = TRUE,  # invierte el enlace log → exp(η) = λ
  re_formula = NA,    # excluye efectos aleatorios
  ndraws      = n_draws
)


dim(lambda_draws)

colnames(lambda_draws)<-paste(newdata$Treatment,newdata$LINE,newdata$SEX)

library(posterior)

shape_draws <- as_draws_df(mw, draws=n_draws)[["shape"]]
shape_draws <- shape_draws[1:n_draws]

TIMEs <- seq(0.1, 7, by = 0.05)  # o el máximo TIEMPO observado en tus datos

ndraws <- nrow(lambda_draws)
ngrupos <- ncol(lambda_draws)
nTIME <- length(TIMEs)

# Crear array vacío
S_arr <- array(NA, dim = c(ndraws, ngrupos, nTIME))

dim(S_arr)

# Calcular S(t) = exp(-(λ * t)^k) para cada combinación
for (g in seq_len(ngrupos)) {
  λ_g <- lambda_draws[, g]  # vector de longitud ndraws
  for (t_idx in seq_len(nTIME)) {
    t <- TIMEs[t_idx]
    S_arr[, g, t_idx] <- exp(- ((t / λ_g)^shape_draws))
  }
}

newdata$group <- seq_len(nrow(newdata))  # asignar índice de grupo
newdata$grupo <- with(newdata, paste(Treatment, LINE, SEX, sep = "_")) 


ndraws <- dim(S_arr)[1]
ngrupos <- dim(S_arr)[2]
ntimes <- dim(S_arr)[3]

df_raw <- tibble(
  draw  = rep(1:ndraws, times = ngrupos * ntimes),
  group = rep(rep(1:ngrupos, each = ndraws), times = ntimes),
  time  = rep(TIMEs, each = ndraws * ngrupos),
  surv  = as.vector(S_arr)
)

df_surv <- df_raw %>%
  left_join(
    newdata[, c("group", "grupo", "Treatment", "LINE", "SEX")],
    by = "group"
  )

df_summary <- df_surv %>%
  group_by(time, grupo, Treatment, LINE, SEX) %>%
  summarise(
    median = median(surv),
    lower  = quantile(surv, 0.05),
    upper  = quantile(surv, 0.95),
    .groups = "drop"
  )


levels(df_summary$Treatment)<- c("Control", "PBO")
levels(df_summary$SEX)<- c("Hembras", "Machos")

windows()
ggplot(df_summary, aes(x = time, y = median, color = Treatment, fill = Treatment)) +
  geom_line(size = 1.2) +
  coord_cartesian(xlim=c(0,max(data$TIME)+0.5))+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3, color = NA,show.legend=F) +
  facet_grid(rows=vars(SEX),cols=vars(LINE)) +
  labs(
    x = "TIME (h)",
    y = "Probabilidad de supervivencia",
    title = ""
  ) +
  guides(color=guide_legend(title="Treatmentamiento"))+
  theme_minimal(base_size = 20)+
  theme(panel.grid=element_blank(), legend.position = "bottom")


windows()
ggplot(df_summary, aes(x = time, y = median, color = LINE, fill = LINE,)) +
  geom_line(size = 1.2) +
  coord_cartesian(xlim=c(0,max(data$TIME)+0.5))+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, color = NA,show.legend=F) +
  facet_grid(rows=vars(SEX),cols=vars(Treatment)) +
  labs(
    x = "TIME (h)",
    y = "Probabilidad de supervivencia",
    title = "Lambda-cialotrina"
  ) +
  guides(color=guide_legend(title="LINE"))+
  theme_minimal(base_size = 20)+
  theme(panel.grid=element_blank(), legend.position = "bottom")

windows()
ggplot(df_summary[df_summary$SEX=="Machos",], aes(x = time, y = median, color = LINE, fill = LINE,)) +
  geom_line(size = 1.2) +
  coord_cartesian(xlim=c(0,max(data$TIME)+0.5))+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, color = NA,show.legend=F) +
  facet_grid(cols=vars(Treatment)) +
  labs(
    x = "TIME (h)",
    y = "Probabilidad de supervivencia",
    title = "Lambda-cialotrina machos"
  ) +
  guides(color=guide_legend(title="LINE"))+
  theme_minimal(base_size = 20)+
  theme(panel.grid=element_blank(), legend.position = "bottom")


windows()
ggplot(df_summary[df_summary$SEX=="Hembras",], aes(x = time, y = median, color = LINE, fill = LINE,)) +
  geom_line(size = 1.2) +
  coord_cartesian(xlim=c(0,max(data$TIME)+0.5))+
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, color = NA,show.legend=F) +
  facet_grid(cols=vars(Treatment)) +
  labs(
    x = "TIME (h)",
    y = "Probabilidad de supervivencia",
    title = "Lambda-cialotrina hembras"
  ) +
  guides(color=guide_legend(title="LINE"))+
  theme_minimal(base_size = 20)+
  theme(panel.grid=element_blank(), legend.position = "bottom")
