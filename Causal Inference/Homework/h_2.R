## homework 2
# setup

library(magrittr)
library(epiDisplay)
library(Hmisc)
library(ggplot2)
library(Matching)
library(overlap)
library(WeightIt)
# import data
setwd('Desktop/USA/GWU/2020FALL/STAT6289 Causal Inference/homework/hw_2')
dat.ful = read.csv('2009Births.csv', header = TRUE)

# b
dat=dat.ful[dat.ful$Smokes==0|dat.ful$Smokes==1,]
dat=dat[complete.cases(dat), ]
dim(dat); n=nrow(dat)
attach(dat)

#”Bmonth”, ”Mage”, ”Hispmom” 
par(mfrow=c(2,1))
hist(dat[dat$Smokes==1,"Meduc"], xlim = c(0,20),main = "'Meduc' in treated group")
hist(dat[dat$Smokes==0,"Meduc"], xlim = c(0,20),main = "'Meduc' in control group")

par(mfrow=c(2,1))
hist(dat[dat$Smokes==1,"Bmonth"],main = "'Bmonth' in treated group")
hist(dat[dat$Smokes==0,"Bmonth"],main = "'Bmonth' in control group")

par(mfrow=c(2,1))
hist(dat[dat$Smokes==1,"Mage"], xlim = c(10,50),main = "'Mage' in treated group")
hist(dat[dat$Smokes==0,"Mage"], xlim = c(10,50),main = "'Mage' in control group")

par(mfrow=c(2,1))
hist(dat[dat$Smokes==1,"Hispmom"]*1,main = "'Hispmom' in treated group")
hist(dat[dat$Smokes==0,"Hispmom"]*1,main = "'Hispmom' in control group")


ggplot(data.frame(dat[dat$Smokes==1,"Hispmom"]), aes(x=dat[dat$Smokes==1,"Hispmom"])) +
  geom_bar()+
  labs(title = '"Hispmom" in treated group')

ggplot(data.frame(dat[dat$Smokes==0,"Hispmom"]), aes(x=dat[dat$Smokes==0,"Hispmom"])) +
  geom_bar()+
  labs(title = '"Hispmom" in control group')


#c
dat1[Gender=='Male','Gender']=1
dat1[Gender == 'Female', 'Gender']=0
dat1[,'Gender'] = as.numeric(dat1[,'Gender'])
dat1[,'Hispmom'] = as.numeric(dat[,'Hispmom'])
dat1[,'Hispdad'] = as.numeric(dat[,'Hispdad'])
attach(dat1)
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
normalized_dif = rep(0,12)
log_ratio_sd = rep(0,12)
pi_treated = rep(0,12)
pi_control = rep(0,12)
tab = data.frame(0, row.names = colnames(dat[,1:12], ))

for (i in 1:12){
  mean_in_treated[i] = mean(dat1[Smokes==1,i])
  mean_in_control[i] = mean(dat1[Smokes==0,i])
  logsd_in_treated[i] = log(sd(dat1[Smokes==1,i]))
  logsd_in_control[i] = log(sd(dat1[Smokes==0,i]))
  normalized_dif[i] = (mean(dat1[Smokes==1,i]) - mean(dat1[Smokes==0,i]))/sqrt((var(dat1[Smokes==1,i])+var(dat1[Smokes==1,i]))/2)
  log_ratio_sd[i] = log(sd(dat1[Smokes==1,i])/sd(dat1[Smokes==0,i])) 
  pi_treated[i] = (sum(dat1[Smokes==1,i]<=min(quantile(dat1[Smokes==1,i], c(0.025, 0.975))))+ sum(dat1[Smokes==1,i]>=max(quantile(dat1[Smokes==1,i], c(0.025, 0.975)))))/nrow(dat1[Smokes==1,])
  pi_control[i] = (sum(dat1[Smokes==0,i]<=min(quantile(dat1[Smokes==0,i], c(0.025, 0.975))))+ sum(dat1[Smokes==0,i]>=max(quantile(dat1[Smokes==0,i], c(0.025, 0.975)))))/nrow(dat1[Smokes==0,])
  }
df = data.frame(mean_in_treated,mean_in_control,
                logsd_in_treated,logsd_in_control,normalized_dif
                ,log_ratio_sd,pi_treated,pi_control)
colnames(df) =c('mean in treated','mean in control',
'log sd. dev in treated','log sd.dev in control',
'normalized difference','log ratio sd. dev.',
'pi in treated','pi in control')
rownames(df) = names(dat)[1:12]
write.table(df,file="c.txt",quote = FALSE, sep = ',')

#d
xbase=c('Marital', 'Meduc')
xnames=names(dat)[1:12] 
xcandidate=setdiff (xnames , xbase)
#linear terms
repeat{
  flag=1 # whether the covariates set is updated 
  for(j in 1:length(xcandidate))
  {
    xnam=c(xbase , xcandidate[j]) 
    model0=glm(as.formula(paste('Smokes~',paste(xbase,collapse='+'))), data=dat1 , family=binomial ( link='logit'))
    model1=glm(as.formula(paste('Smokes~',paste(xnam,collapse= "+"))), data=dat1 , family=binomial ( link='logit'))
    if (lrtest (model0 , model1)$p.value <0.05){ 
      xbase=xnam
      flag=0
    } }
  xcandidate=setdiff(xcandidate ,xbase)
  if(flag) break }
xbase # includes all selected linear terms

#TotPreg, Visits, Hispmon,Hispdad Feduc

#e

model_1=glm(Smokes~(Marital+Meduc+TotPreg+Visits+Hispdad+Hispmom+Feduc)^2+I(Marital^2)+I(Meduc^2)+I(TotPreg^2)+I(Visits^2)+I(Hispdad^2)+I(Hispmom^2)+I(Feduc^2),data = dat1,family=binomial(link='logit'))
pv = coef(summary(model_1))[,4]
pv = data.frame(pv[pv<0.05])
a = row.names(pv)
a
#"Marital"       "Visits"        "I(TotPreg^2)"  "Meduc:Hispmom" "Meduc:Feduc"  

#f
model_new = glm(Smokes~ Marital + Meduc + TotPreg + Hispdad + Hispmom + Feduc + Visits + Marital:Marital+ I(TotPreg^2)+ Meduc:Hispmom + Meduc:Feduc, data = dat1,
                family=binomial(link='logit'))
summary(model_new)
dat2= dat1
dat2 = cbind(dat2, ps = predict(model_new, type="response", newdata=dat2))
dat2 = cbind(dat2, lps = predict(model_new, newdata=dat2))
hist(dat2[dat2$Smokes==1,'lps'],main = "propensity score in treated",xlab = 'linearlized propensity score')
hist(dat2[dat2$Smokes==0,'lps'],main = "propensity score in control",xlab = 'linearlized propensity score')

par(mfrow=c(1,1))
histbackback(split(dat2$lps,	dat2$Smokes),	main= "Propensity	score	before	matching",	xlab=c("control",	"treatment"))

#g
#create index collumn
min(dat2[dat2$Smokes==1,'lps'])
max(dat2[dat2$Smokes==0,'lps'])

o = data.frame(subset(dat2, Smokes==1 & lps> 0.1045609))
p = data.frame(subset(dat2, Smokes==0 & lps< (-5.007004)))
o_row = rownames(o)
p_row = rownames(p)
del_index = as.numeric(append(o_row,p_row))
trimmed.dat = dat2[-del_index,]
detach(dat2)
attach(trimmed.dat)
nrow(trimmed.dat[Smokes==1,])
#189
nrow(trimmed.dat[Smokes==0,])
#1347
par(mfrow=c(2,1))
hist(trimmed.dat[trimmed.dat$Smokes==1,'lps'],main = "propensity score in treated after trimming(treated)",xlab = 'linearlized propensity score')
hist(trimmed.dat[trimmed.dat$Smokes==0,'lps'],main = "propensity score in control after trimming(control)",xlab = 'linearlized propensity score')

#h
J = 1
K = 0
trimmed.dat$block = 1
while(J>K){
  for (i in (K+1):J){
    dat.c = trimmed.dat[trimmed.dat$block==i,]
    lps.c = dat.c$lps
    dat_a = dat.c[lps.c<median(lps.c),]
    dat_b = dat.c[lps.c>=median(lps.c),]
    tt = t.test(dat_a$lps, dat_b$lps)
    if((tt$statistic < 0.05)){
      ntl = sum(dat_a$Smokes==1);ncl=sum(dat_a$Smokes==0);nl=ntl+ncl
      ntu = sum(dat_b$Smokes==1);ncu=sum(dat_b$Smokes==0);nu=ntu+ncu
      if((min(ntl,ncl,ntu,ncu)>=3) & (min(nl,nu)>=p+2)){
        dat.c$block[lps.c<median(lps.c)] = J+1
        dat.c$block[lps.c>=median(lps.c)] = J+2
        trimmed.dat$block[trimmed.dat$block==i]=dat.c$block
        trimmed.dat$block = as.numeric(as.factor(trimmed.dat$block))
        J = J+1
      }else{K=K+1}
    }else{K=K+1}
  }
}
J
# J = 7

#matching
matched.dat = trimmed.dat[order(lps, decreasing = TRUE),]

#i
attach(trimmed.dat)
norm_dif_trimmed = rep(0,12)
dif = data.frame(matrix(0,nrow = 12,ncol = 7))
n = nrow(trimmed.dat)
nj = rep(0,7)
norm_dif_block = rep(0,12)
for (i in 1:12) {
  norm_dif_trimmed[i] = (mean(trimmed.dat[Smokes==1,i]) - mean(trimmed.dat[Smokes==0,i]))/sqrt((var(trimmed.dat[Smokes==1,i])+var(trimmed.dat[Smokes==0,i]))/2)
  for (j in 1:7){
    dif[i,j] = (mean(trimmed.dat[Smokes==1&block==j,i])-mean(trimmed.dat[Smokes==0&block==j,i]))/sqrt((var(trimmed.dat[Smokes==1&block==j,i]) + var(trimmed.dat[Smokes==0&block==j,i]))/2)
    nj[j] = nrow(trimmed.dat[block==j,])
    if (dif[i,j] =='-Inf' | dif[i,j] == 'NaN'){
      dif[i,j] = 0
    }
    norm_dif_block[i] = norm_dif_block[i] + (dif[i,j] * (nj[j]/n))
  }
}

norm_dif = cbind(normalized_dif,norm_dif_trimmed, norm_dif_block)
rownames(norm_dif) = names(dat)[1:12]
write.table(norm_dif,file="i.txt",quote = FALSE, sep = ',')


#j
z = data.frame(matrix(0,nrow = 7,ncol = 12))
for (j in 1:7) {
  for (k in 1:12) {
    xt = mean(trimmed.dat[Smokes==1&block==j,k])
    xc = mean(trimmed.dat[Smokes==0&block==j,k])
    nc = nrow(trimmed.dat[Smokes==0 & block==j,])
    nt = nrow(trimmed.dat[Smokes==1 & block==j,])
    deno = sqrt(var(trimmed.dat[block==j,k])*(1/nc + 1/nt))
    z[j,k] = (xt-xc)/deno
  }
  
}

par(mfrow=c(1,1))
z1 = unlist(z)
qqnorm(z1)
qqline(z1,col=2)

#k

tau_ht = (1/n) *sum(((Smokes-ps)*BirthWeight)/(ps*(1-ps)))
lambda_ht = abs((Smokes - ps)/(ps*(1-ps)))
max(lambda_ht)
min(lambda_ht)
sd(lambda_ht)
ht = c(max(lambda_ht),min(lambda_ht),sd(lambda_ht))


lambda_strat = rep(0, nrow(trimmed.dat))
for (i in 1:n) {
  for (j in 1:7) 
    {
    if (trimmed.dat[i,'block']==j & trimmed.dat[i,'Smokes']==0) {
      lambda_strat[i] = lambda_strat[i] + nj[j]/(nrow(trimmed.dat[Smokes==0&block==j,]))
    }
    else if (trimmed.dat[i,'block']==j & trimmed.dat[i,'Smokes']==1)
      {
      lambda_strat[i] = lambda_strat[i] + nj[j]/(nrow(trimmed.dat[Smokes==1&block==j,]))
    }
  }
  }
max(lambda_strat)
min(lambda_strat)
sd(lambda_strat)
strat = c(max(lambda_strat),min(lambda_strat),sd(lambda_strat))
tau_strat = 1/n * (sum(Smokes*BirthWeight*lambda_strat) - sum((1-Smokes)*BirthWeight*lambda_strat))

ht_vs_srtat = data.frame(ht,strat,row.names = c('Maximum','Minimum','Standard Deviation'))

#i
tau_ht
tau_strat

fmodel = glm(BirthWeight~ I(tau_strat*Smokes))
summary(fmodel)

#m
variance = 0
for (i in n) {
  if (trimmed.dat[i,'Smokes']==1){
  variance = variance +(1/(n^2)) * (lambda_ht^2)[i] * (var(trimmed.dat[Smokes==1,i]))
  }else if(trimmed.dat[i,'Smokes']==0){
  variance =  variance + (1/(n^2)) * (lambda_ht^2)[i] * (var(trimmed.dat[Smokes==0,i]))
  }
}







