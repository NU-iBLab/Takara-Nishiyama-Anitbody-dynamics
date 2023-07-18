rm(list=ls(all=TRUE))

CURRENT_WORKING_DIR <- dirname(rstudioapi::getActiveDocumentContext()$path)
CURRENT_WORKING_DIR
setwd(CURRENT_WORKING_DIR)

###################################################
## Install package
###################################################
{
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
  library(deSolve)
  library(FME)
  library(patchwork)
  library(fastDummies)
  library(ropls)
  library(reticulate)
}
###################################################
## function
###################################################

## fatal ##
system("R CMD SHLIB dec_model.c")
dyn.load( paste("dec_model",.Platform$dynlib.ext,sep="") )

## ODE model ##
ODEs <- function(pars) {
  
  #initial value is squared parameters
  rhs <- c(M=M0)
  times <- seq(Tmin,Tmax,step_size)
  
  ## C compiled version ##
  pars <- c(pars)
  out <- ode(y=rhs,parms=pars,times=times,func="derivs",initfunc="initparms",nout=1,outnames=c(""),dllname="dec_model",method="rk4")
  as.data.frame(out)
}

###################################################
## Data setting
###################################################

df_par <- read.csv("../monolix/estimatedIndividualParameters.txt")[,c(1,19:25)]
df_obs <- read.csv("../monolix/sample_observation.csv")
par.set <- left_join(df_par,unique(df_obs[,c(1,5,7)]),by="id")
colnames(par.set)[2:10] <- c("M0","M33","mu1","mu2","m","K","delay3","vac3","group")
exp.plot <- df_obs[,c(1,2,3,5,7)] %>% 
  mutate(categorical=as.factor(categorical),
         time_b=time-vac3) %>% 
  select("id","time_b","obs","categorical")
colnames(exp.plot)[c(2,4)] <- c("time","group")

###################################################
## Analysis
###################################################
id <- unique(par.set$id)

Tmin <- 0
Tmax <- 1500.0
step_size <- 1.0
stime <- seq(Tmin,Tmax,step_size)

simu.df <- c()
peak.list <- rep(NA,length(id))
duration <- rep(NA,length(id))
uprate <- rep(NA,length(id))
A0 <- rep(NA,length(id))
for (i in 1:length(id)){
  pars <- par.set[par.set$id==id[i],3:9] %>% as.numeric()
  M0 <- par.set$M0[i]
  
  fitted <- ODEs(pars)
  
  peak.list[i] <- max(fitted[fitted$time>=pars[7],]["M"])
  a <- fitted[fitted$time>=pars[7] & fitted$M>80,]
  duration[i]<-a$time[length(a$time)]-a$time[1]
  A0[i] <- fitted$M[pars[7]+1]
  amin <- (fitted[fitted$time>=(pars[7]+pars[6]-1),]["M"][1,])
  uprate[i] <- peak.list[i]/amin
  
  kari.df <- data.frame(id=rep(id[i],nrow(fitted)),
                        time=fitted$time-par.set$vac3[i],
                        value=fitted$M)
  
  simu.df <- rbind(simu.df,kari.df)
}
static.par <- cbind(par.set[,1:8],A0=A0,peak=peak.list,duration=duration,
               uprate=uprate,group=par.set$group)

### PLS-DA ###
simu.df.l <- simu.df[simu.df$time>=-250 & simu.df$time<=430,]
ClusNum <- data.frame(id=par.set$id,group=as.numeric(par.set$group)+1)
dummy2 <- dummy_cols(ClusNum$group)

df <- reshape(simu.df.l, idvar = "id", timevar = "time", direction = "wide")

vlRF1.plsda <- opls(df[,2:682], dummy2$.data_1)
vlRF2.plsda <- opls(df[,2:682], dummy2$.data_2)

mean <- colMeans(df[,2:682])
sd <- apply(df[,2:682], 2, sd)

R1 <- mean + sd*(mean(vlRF1.plsda@scoreMN[dummy2$.data_1==1,1]))*vlRF1.plsda@loadingMN[,1]
R2 <- mean + sd*(mean(vlRF2.plsda@scoreMN[dummy2$.data_2==1,1]))*vlRF2.plsda@loadingMN[,1]

R1_1 <- mean+sd*(mean(vlRF1.plsda@scoreMN[dummy2$.data_1==1,1])+sd(vlRF1.plsda@scoreMN[dummy2$.data_1==1,1]))*vlRF1.plsda@loadingMN[,1]
R1_2 <- mean+sd*(mean(vlRF1.plsda@scoreMN[dummy2$.data_1==1,1])-sd(vlRF1.plsda@scoreMN[dummy2$.data_1==1,1]))*vlRF1.plsda@loadingMN[,1]

R2_1 <- mean+sd*(mean(vlRF2.plsda@scoreMN[dummy2$.data_2==1,1])+sd(vlRF2.plsda@scoreMN[dummy2$.data_2==1,1]))*vlRF2.plsda@loadingMN[,1]
R2_2 <- mean+sd*(mean(vlRF2.plsda@scoreMN[dummy2$.data_2==1,1])-sd(vlRF2.plsda@scoreMN[dummy2$.data_2==1,1]))*vlRF2.plsda@loadingMN[,1]

limit <- -0.5

df_R1 <- data.frame(R1) 
df_R2 <- data.frame(R2) 

df_R1$time <- seq(-250,430,1)
df_R1$lower <- R1_1
df_R1$upper <- R1_2

df_R2$time <- seq(-250,430,1)
df_R2$lower <- R2_1
df_R2$upper <- R2_2

df_R1$lower[df_R1$lower < limit] <- limit
df_R1$upper[df_R1$upper < limit] <- limit
df_R2$lower[df_R2$lower < limit] <- limit
df_R2$upper[df_R2$upper < limit] <- limit


### Statistic Analysis###
static.par.0 <- static.par[static.par$group=="0",]
static.par.1 <- static.par[static.par$group=="1",]

parname <- colnames(static.par)
p.res <- rep(NA,length(id))
for (i in 2:12) {
  res <- wilcox.test(static.par.0[,i],static.par.1[,i])
  if(res$p.value>0.05){
    icon <- "NS."
  }else if(res$p.value>0.01){
    icon <- "*"
  }else if(res$p.value>0.001){
    icon <- "**"
  }else if(res$p.value>0.0001){
    icon <- "***"
  }else{
    icon <- "****"
  }
  msd_n <- paste0(round(mean(static.par.0[,i]),4),"(",round(sd(static.par.0[,i]),4),")",sep=",")
  msd_r <- paste0(round(mean(static.par.1[,i]),4),"(",round(sd(static.par.1[,i]),4),")",sep=",")
  d.res <- data.frame(par=parname[i],
                      naive=msd_n,
                      recovered=msd_r,
                      p.val=res$p.value,pcol=icon)
  p.res <- rbind(p.res,d.res)
}
