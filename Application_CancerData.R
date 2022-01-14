# NOTATION
# Let Yijl be the random variable corresponding to the digital RNA expression (DRE) measure
# for observation i (i = 1, . . . , n with n = 714)
# of condition j (j = 0, . . . , d with d = 3)
# in biological replicate l (l = 1, . . . , rj)
# with yijl being the corresponding observed value of Yijl .

# j = 4 conditions : S, ADC, ASC, SCC
# 6  = r1	replicate ADS		(cond1)
# 2  = r2	replicate ASC		(cond2)
# 21 = r3	replicate SCC		(cond3)
# 29 = r0	replicate S 		(cond4)

# We define an observation as a biological entity, such as a gene or RNA,and a variable as a
# replicate in a given condition.
# Let q = sum_j rj = 29 + 6 + 2 + 21 = 58 be the total number of variables (all replicates 
# in all conditions) in the data, such that bold{y}= (y_ijl) 714x58 matrix of the DRE for all
# observations and variables, and yi is the q-dimensional vector of DRE for all variables of
# observation i with y_i belonging to N^q
# We use dot notation to indicate summations in various directions, e.g.,
# y_.jl = sum_i y_ijl (sum by sRNA - col)
# y_i.. = sum_j sum_l y_ijl (sum by condition and replicate - row)
# Let s_jl represent the normalized library size for replicate l of condition j.

# LOAD PACKAGES
library(MASS)
library(deepMOU)
library(MM)
library(MGLM)
source("FS_MorelNagaraj.R")
source("EM_mixtures DM.R")

# LOAD DATA
# From "Ultra-high throughput sequencing-based small RNA discovery and discrete statistical
# biomarker analysis in a collection of cervical tumours and matched controls" Witten e al. (2010)
setwd("~/Dropbox/Tesi Noemi/per analisi dati reali")
x <- read.csv("GSE20592.csv", dec=",", sep=";", header=TRUE)

#varNames <- c("N1", "T1", "N2", "T2", "N3", "T3", "N4", "T4", "N5", "T5", "N6", "T6", "N7", "T7", "N8", "T8", "N9", "T9", "N10", "T10", "N11", "T11", "N12", "T12", "N13", "T13", "N14", "T14", "N15", "T15", "N16", "T16", "N17", "T17", "N18", "T18", "N19", "T19", "N20", "T20", "N21", "T21", "N22", "T22", "N23", "T23", "N24", "T24", "N25", "T25", "N26", "T26", "N27", "T27", "N28", "T28", "N29", "T29")
#colnames(x) <- varNames
x[is.na(x)] <- 0  #Sostitution of NA with 0
x = x[rowSums(x) > 0,]
x = x[,colSums(x) > 0]
rownames(x) <- paste("sRNA", 1:nrow(x), sep="_")
head(x)


# The data consist of counts for 714 sRNA 
# 714 obs. of  58 variables
summary(x)
max(x)  # The counts range in size from 0 to 476438
table(x[x==0])	# 19733 null obs.

# CREATE 2 DATASET: Normal and Tumor ---> DEVO TRASPORRE IL DATASET????
# ADC, adenocarcinoma = 6  replicates -> col. 2, 4, 6, 8, 10, 12 
# ASC, adenosquamous cell carcinoma	= 2  replicates -> col. 14, 16
# SCC, squamous cell carcinoma = 21 replicates -> col. 18, 20, 22, etc..
# Normal, without carcinoma = 29 replicates -> col. 1, 3, 5, 7, etc..
ASC = cbind(x$G871T, x$G701T)
ADC = cbind(x$G547T, x$G659T, x$G691T, x$G696T, x$G761T, x$G220T)
SCC = cbind(x$G576T, x$G652T, x$G699T, x$G850T, x$G013T, x$G603T, 
            x$G026T, x$G575T, x$G613T, x$G702T, x$G623T, x$G645T, 
            x$G529T, x$G648T, x$G601T, x$G727T, x$G001T, x$G243T, 
            x$G531T, x$G612T, x$G428T)

Tumor = cbind(ADC, SCC, ASC)
Normal = cbind(x$G547N, x$G659N, x$G691N, x$G696N, x$G761N, x$G220N, x$G576N, x$G652N, x$G699N, 
             x$G850N, x$G013N, x$G603N, x$G026N, x$G575N, x$G613N, x$G702N, x$G623N, x$G645N, 
             x$G529N, x$G648N, x$G601N, x$G727N, x$G001N, x$G243N, x$G531N, x$G612N, x$G428N, 
             x$G871N, x$G701N)
Tumor <- t(Tumor)
Normal <- t(Normal)



############################### ESTIMATES MODELS #####################################
data = Normal # Or Tumor
data = rbind(Normal,Tumor)
table(is.na(data))
data = data[rowSums(data) > 0,]
data = data[,colSums(data) > 0]
maxdata=max(data)
data=data/maxdata

# MULTINOMIAL
fit0 = MGLMfit(data, dist="MN")
bic_mn = fit0@BIC


# DIRICHLET-MULTINOMIAL
#fit = MGLMfit(data, dist="DM") -- Con questa non fitta, è lentissimo
#sum_alpha = sum(fit@estimate)
#theta = fit@estimate/sum_alpha
#rho = 1/sqrt(1 + sum_alpha) 
#bic_dm = fit@BIC
fit = dir_mult_GD(data, 1)
sum_alpha = sum(fit$Theta)
theta = fit$Theta/sum_alpha
rho = 1/sqrt(1 + sum_alpha)
llik_dm = fit$likelihood[length(fit$likelihood)]
n = nrow(data)
p = ncol(data)
bic_dm = -2*llik_dm + p*log(n)

# NEERCHAL&MOREL ## lentissimo
fit2 = FS_mix_mult(data, rho = rho, pi_greco = theta)
rho2 = fit2$rho
theta2 = fit2$pi_greco
llik_rcm = fit2$lik_theta[length(fit2$lik_theta)]
n = nrow(data)
p = ncol(data)
bic_rcm = -2*llik_rcm + p*log(n)
aic_rcm= -2*llik_rcm + 2*p

# NEGATIVE MULTINOMIAL # DA PROBLEMI
fit3 = MGLMfit(data, dist="NegMN")
pp = ncol(data)
prob = fit3@estimate[1:pp]
prob0 = 1 - sum(prob)
phi = fit3@estimate[pp+1]
bic_nm = fit3@BIC

# GENERALIZED DIRICHLET MULTINOMIAL # NON VA
B = 0
Bi = 0
fit4 = MGLMfit(data, dist="GDM")
alpha = fit4@estimate[1:(pp-1)]
beta2 = fit4@estimate[pp:length(fit4@estimate)]
bic_gdm = fit4@BIC
A = alpha/(alpha + beta2)
Ai = (alpha + 1)/(alpha + beta2 + 1)
B = cumprod(beta2/(beta2 + alpha))
Bi = cumprod((beta2 + 1)/(beta2 + alpha + 1))
B = c(1, B[-pp])
Bi = c(1, Bi[-pp])

# DEEP DIRICHLET-MULTINOMIAL
# K=2
fit5_2 = em.mixt.DM(data, k=2, it=200, eps=1e-5, seme=7, KK=10, qq=10)
bic_ddm_2 = fit5_2$bic

# K=3
fit5_3 = em.mixt.DM(data, k=3, it=200, eps=1e-5, seme=7, KK=10, qq=10)
bic_ddm_3 = fit5_3$bic

# K=4
fit5_4 = em.mixt.DM(data, k=4, it=200, eps=1e-5, seme=7, KK=10, qq=10)
bic_ddm_4 = fit5_4$bic

# K=20
fit5_20 = em.mixt.DM(data, k=20, it=200, eps=1e-5, seme=7, KK=10, qq=10)
bic_ddm_20 = fit5_20$bic

# COMPUTED VARIANCE
var_comp = apply(data, 2, var)
var_comp = var_comp*maxdata^2 ## ritorno alle varianze dei dati orginali

# ESTIMATED PROBABILITY
prob = colSums(data)/sum(data)

# ESTIMATED VARIANCES
m = mean(rowSums(data))
var_mult = m%*%t(prob*(1-prob))*maxdata^2
var_dmult = ((1+rho^2*(m-1))*m)%*%t(theta*(1-theta))*maxdata^2
var_rcm = ((1+rho2^2*(m-1))*m)%*%t(theta2*(1-theta2))*maxdata^2
var_deepmult2 = mean(diag(fit5_2$varianza))*maxdata^2
var_deepmult3 = mean(diag(fit5_3$varianza))*maxdata^2
var_deepmult4 = mean(diag(fit5_4$varianza))*maxdata^2
var_deepmult20 = mean(diag(fit5_20$varianza))*maxdata^2
#var_negmult = phi*prob/prob0*(1 + prob/prob0)
#var_gdm = m*A*B*((m-1)*Bi*Ai - m*B*A + 1)


# EUCLIDEAN DISTANCES (divided by 1000000)
sqrt(mean((var_comp-var_mult)^2))/1000000
sqrt(mean((var_comp-var_dmult)^2))/1000000
sqrt(mean((var_comp-var_rcm)^2))/1000000
#100*sqrt(sum((var_comp-var_negmult)^2))
#100*sqrt(sum((var_comp-var_gdm)^2, na.rm = TRUE))
sqrt(mean((var_comp-var_deepmult2)^2))/1000000
sqrt(mean((var_comp-var_deepmult3)^2))/1000000
sqrt(mean((var_comp-var_deepmult4)^2))/1000000
sqrt(mean((var_comp-var_deepmult20)^2))/1000000

save.image("analisi.RData")


## il bic non ci sceglie quindi non metterei l'enfasi a questo
# BIC
bic_mn
bic_dm
bic_rcm
bic_nm
min(bic_mn, bic_dm, bic_rcm, bic_nm)
bic_gdm
bic_ddm_2
bic_ddm_3
bic_ddm_4
bic_ddm_20



alpha=fit5_3$alpha

























# NEERCHAL&MOREL
fit2 = try(FS_mix_mult(x.zero, rho=rho, pi_greco = theta))
if (!is.character(fit2)) {
  rho2 = fit2$rho
  theta2 = fit2$pi_greco
  llik_rcm = fit2$lik_theta[length(fit2$lik_theta)]
  bic_rcm[i,j] = -2*llik_rcm + p*log(n)
  aic_rcm[i,j] = -2*llik_rcm + 2*p}

# NEGATIVE MULTINOMIAL
fit3 = try(MGLMfit(x.zero, dist="NegMN"))
pp=ncol(x.zero)
if (!is.character(fit3)) {
  prob = fit3@estimate[1:pp]
  prob0 = 1 - sum(prob)
  phi = fit3@estimate[pp+1]
  bic_nm[i,j] = fit3@BIC
  aic_nm[i,j] = fit3@AIC}

# GENERALIZED DIRICHLET MULTINOMIAL
B=0
Bi=0
fit4 = try(MGLMfit(x.zero, dist="GDM"))
if (!is.character(fit4)) {
  alpha = fit4@estimate[1:(pp-1)]
  beta2 = fit4@estimate[pp:length(fit4@estimate)]
  bic_gdm[i,j] = fit4@BIC
  aic_gdm[i,j] = fit4@AIC
  A = alpha/(alpha + beta2)
  Ai = (alpha + 1)/(alpha + beta2 + 1)
  B=cumprod(beta2/(beta2+alpha))
  Bi=cumprod((beta2+1)/(beta2+alpha+1))
  B=c(1,B[-pp])
  Bi=c(1,Bi[-pp])}

# MISTURE 3
fit5_2 = em.mixt.DM(x.zero, k=2, it=200, eps=1e-5, seme=7, KK=10, qq=10)
bic_ddm2[i,j] = fit5_2$bic
aic_ddm2[i,j] = fit5_2$aic

# MISTURE 5
fit5_3 = em.mixt.DM(x.zero, k=3, it=200, eps=1e-5, seme=7, KK=10, qq=10)
bic_ddm3[i,j] = fit5_3$bic
aic_ddm3[i,j] = fit5_3$aic

# MISTURE 6
fit5_4 = em.mixt.DM(x.zero, k=4, it=200, eps=1e-5, seme=7, KK=10, qq=10)
bic_ddm4[i,j] = fit5_4$bic
aic_ddm4[i,j] = fit5_4$aic

# MISTURE 20
fit5_20 = em.mixt.DM(x.zero, k=20, it=200, eps=1e-5, seme=7, KK=10, qq=10)
bic_ddm20[i,j] = fit5_20$bic
aic_ddm20[i,j] = fit5_20$aic































