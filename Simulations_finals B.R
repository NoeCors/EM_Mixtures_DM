###################### MULTINOMIAL #######################
# The 50 samples drawn are in COLUMNS. Rows are the CATEGORIES. 
# Size is the sum of all the outcomes/experiments in/for each category.
# The function rmultinom gives as result a Kxn matrix.
#library(devtools)
#install_github("andrewraim/COMMultReg")
#library(COMMultReg)
#source("ZIGDM_EM.R")

# Simulazione1 = Random Zeros
# Simulazione2 = Second way of inserting zeros
# Simulazione3 = Proof of convergence of the variance of the MIXTURE

library(MASS)
library(deepMOU)
library(MM)
library(MGLM)
source("FS_MorelNagaraj.R")
source("EM_mixtures DM.R")


  p = 10 
  n = 200 
  nrep = 10 
  qzeros = seq(0, 1, length.out = 11)[-11]
  var_orig_comp = var_zero_comp = matrix(0,nrep,length(qzeros))  ## prendo varianza media
  var_mult_x = var_mult_a = var_dmult_a = var_nm_a = var_deepmult_a2 = var_deepmult_a3 = var_deepmult_a4 = matrix(NA, nrep, length(qzeros))
  var_deepmult_a20 = var_negmult_a = var_gdm_a = matrix(NA, nrep, length(qzeros))
  bic_dm = bic_mn = bic_nm = bic_gdm = bic_ddm20 = bic_rcm = bic_ddm2 = bic_ddm3 = bic_ddm4 = matrix(NA, nrep, length(qzeros))
  aic_dm = aic_mn = aic_nm = aic_gdm = aic_ddm20 = aic_rcm = aic_ddm2 = aic_ddm3 = aic_ddm4 = matrix(NA, nrep, length(qzeros))
  
  seed = 0
  
  for (i in 1:nrep) for (j in 1:10) { 
    
    print(paste(i,j,sep=" "))
    set.seed(i)
    prob_true = runif(p)
    prob_true = prob_true/sum(prob_true)
    size = 5*p
    x = t(rmultinom(n, size, prob = prob_true))
    xmedia=round(mean(x))
    
    # ADDING ZEROS
    x.zero = x
    x.zero = ifelse(x < (j+2), 0, x) # Simulazione2
    index.row=which(rowSums(x.zero) == 0)
    index.column=which(colSums(x.zero) == 0)
    if (length(index.row)>0) x.zero[index.row,sample(1:p,1)]=xmedia
    if (length(index.column)>0) x.zero[sample(1:n,1),index.column]=xmedia
    
    
    # MULTINOMIAL
    fit0 = try(MGLMfit(x.zero, dist="MN"))
    if (!is.character(fit0)) {
    bic_mn[i,j] = fit0@BIC
    aic_mn[i,j] = fit0@AIC}
    
    # DIRICHLET-MULTINOMIAL
     fit = try(MGLMfit(x.zero, dist="DM"))
     if (!is.character(fit)) {
     sum_alpha = sum(fit@estimate)
     theta = fit@estimate/sum_alpha
     rho = 1/sqrt(1 + sum_alpha) 
     bic_dm[i,j] = fit@BIC
     aic_dm[i,j] = fit@AIC}
    
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

    # COMPUTED VARIANCES
    var_orig_comp[i,j] = mean(apply(x, 2, var))
    var_zero_comp[i,j] = mean(apply(x.zero, 2, var))
    
    # ESTIMATED PROBABILITIES
    prob_s = colSums(x)/sum(x)
    prob_s_zeros = colSums(x.zero)/sum(x.zero)
   
    # ESTIMATED VARIANCES
    m = mean(size)
    m.zero = mean(rowSums(x.zero))
    var_mult_x[i,j] = mean(m*prob_s*(1-prob_s))
    var_mult_a[i,j] = mean(m.zero%*%t(prob_s_zeros*(1-prob_s_zeros)),na.rm=TRUE)
    var_dmult_a[i,j] = mean(((1+rho^2*(m.zero-1))*m.zero)%*%t(theta*(1-theta)),na.rm=TRUE)  
    var_nm_a[i,j] = mean(((1+rho2^2*(m.zero-1))*m.zero)%*%t(theta2*(1-theta2)),na.rm=TRUE)
    var_deepmult_a2[i,j] = mean(diag(fit5_2$varianza),na.rm=TRUE)
    var_deepmult_a3[i,j] = mean(diag(fit5_3$varianza),na.rm=TRUE)
    var_deepmult_a4[i,j] = mean(diag(fit5_4$varianza),na.rm=TRUE)
    var_deepmult_a20[i,j] = mean(diag(fit5_20$varianza),na.rm=TRUE)
    var_negmult_a[i,j] = mean(phi*prob/prob0*(1 + prob/prob0),na.rm=TRUE)
    if (!is.character(fit4)) var_gdm_a[i,j] = mean(m.zero*A*B*((m.zero-1)*Bi*Ai - m.zero*B*A + 1),na.rm=TRUE)
    }
  
  save.image("SimulazioneB_AIC_BIC.RData") 
                                   
  