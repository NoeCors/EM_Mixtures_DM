
em.mixt.DM=function(y,k,it=200,eps=1e-5,seme=7,KK=30,qq=50)
{
  # y sono i dati
  # k il numero di livelli della variabile latente delle componenti
  # it il numero delle iterazioni
  
  set.seed(seme)
  numobs=nrow(y)
  p=ncol(y)
  
  #####################################################
  ## inizializzazione 
    
    
    ## inizializzazione dei pesi 
    pi.greco=rep(1/k,k)
    
    ## inizializzazione di alpha
    alpha=matrix(runif(p*k,-10,10),p,k)  
    for (i in 1:k) alpha[,i]=norm_alpha(alpha[,i])
    
   ## inizializzazione di Beta
    Beta=colMeans(y)      
    Beta=ifelse(Beta==0,exp(-500),Beta)      
    
  
    f.z.y<-f.y.z<-matrix(0,numobs,k)
  
  
  likelihood<-NULL
  ratio=1000
  lik=-10^10
   h=0
  
  while ((h < it) & (ratio > eps )) {
    h=h+1
    #print(h)
    ######## E-STEP 
  # Stima di f(y|z) e f(z|y)
    
      for (i in 1:k) {val=calc.f.y.z(y,alpha[,i],Beta)
                      if (min(val) < (-740)) val=val -740 - min(val) 
                      f.y.z[,i]<-exp(val)
                      f.z.y[,i]<-(f.y.z[,i]*pi.greco[i])}  
    

    f.y<-f.y.z%*%pi.greco
    # Stima di f(z|y): normalizzo per il denominatore 
    f.z.y<-f.z.y/matrix(rowSums(f.z.y),numobs,k,byrow=FALSE)
    
    
    ######## M-STEP 
    
    # Stima pigreco
    pi.greco<-apply(f.z.y,2,mean)

    
    # Stima di Beta
        
      Beta<-nlminb(Beta,beta.f,gradient=beta.f.grad,lower=exp(-700),upper=100,numobs=numobs,k=k,y=y,alpha=alpha,f.z.y=f.z.y, control = list(iter.max = KK))$par
      
     #  for (h in 1:KK) {
     #    out=beta.f.grad(Beta,numobs,k,y,alpha,f.z.y)
     #  Beta=Beta-c(t(out$grad)%*%ginv(out$hess))
     #  }
     #  ## c'è un problema di scala 
     # Beta=Beta/max(Beta)
      
                    
    # Stima di alpha
        
      for (i in 1:k){ 
        alpha[,i]=nlminb(alpha[,i],alpha.f,gradient=alpha.f.grad,lower=-0.9999,upper=0.9999,y=y,Beta=Beta,f.z.y=f.z.y[,i], control = list(iter.max = KK))$par
          }
        
        
      temp=sum(log(f.y))
      
      likelihood=c(likelihood,temp) 
      ratio<-(temp-lik)/abs(lik)
      if (h < qq) ratio=2*eps  ## questo serve a garantire che vengano fatte almeno qq iterazioni
      lik<-temp
      
      if (is.na(lik)) ratio<-eps/2
      
      }
  
  
  hh=(k-1)+p*k+p
  
  lik.ultima=temp
  bic= -2*lik.ultima+hh*log(numobs)
  aic= -2*lik.ultima+2*hh
  
  
  # classificazione 
  cl.k<-apply(f.z.y,1,which.max)
   
  # calcolo della varianza
  
  m=mean(rowSums(y))
  theta=matrix(Beta,p,k)*(1+alpha)
  theta0=colSums(theta)
  lambda=theta/matrix(theta0,p,k,byrow=TRUE)
  rho.quadro=1/(1+theta0)
  within=matrix(0,p,p)
  between=matrix(0,p,p)
  media=m*rowSums(matrix(pi.greco,p,k,byrow=TRUE)*lambda)
  for (i in 1:k) {within=within + pi.greco[i]*m*(diag(lambda[,i])-lambda[,i]%*%t(lambda[,i]))*(1+rho.quadro[i]*(m-1))
                  between=between + pi.greco[i]*(m*lambda[,i]-media)%*%t(m*lambda[,i]-media)
                  }
  varianza=within+between
  
  output=list(likelihood=likelihood,cl=cl.k,pi.greco=pi.greco,bic=bic,aic=aic,Beta=Beta,alpha=alpha,varianza=varianza)
  invisible(output)
}



## Funzioni per stima di Beta #####

 # beta: funz. da massimizzare ####
 beta.f<-function(Beta,numobs,k,y,alpha,f.z.y){
   p=ncol(y)
   temp1<-matrix(NA,numobs,k)
   for (i in 1:k) temp1[,i]=calc.f.y.z(y,alpha[,i],Beta)
   return(-sum(temp1*f.z.y))
 }  



 # calc.f.y.z #### 
calc.f.y.z<-function(y,alpha,Beta){
  numobs=nrow(y)
  p=ncol(y)
  alpha.star=Beta*(1+alpha)
  alpha.star=ifelse(alpha.star==0,0.0001,alpha.star)
  alpha.star.val=rep(t(Beta)%*%(1+alpha),numobs)
  f.y.z<-lgamma(alpha.star.val)-lgamma(rowSums(y)+alpha.star.val)+rowSums(lgamma(y+matrix(alpha.star,numobs,p,byrow=T)))-rowSums(lgamma(matrix(alpha.star,numobs,p,byrow=T)))
  return(f.y.z)
}

# beta: funz. per il gradiente #####
beta.f.grad<-function(Beta,numobs,k,y,alpha,f.z.y){
  p=ncol(y)
  temp1<-matrix(0,numobs,p)
  for (i in 1:k) temp1=temp1+matrix(f.z.y[,i],numobs,p)*calc2grad.b(y,alpha[,i],Beta)
  #return(list(grad=-apply(temp1,2,sum),hess=t(temp1)%*%(temp1)))
  return(-apply(temp1,2,sum))
} 
 

 
# # beta: calcolo gradiente #####
 calc2grad.b<-function(y,alpha,Beta){
   numobs=nrow(y)
   p=ncol(y)
   alpha.star.v=Beta*(1+alpha)
   alpha.star=c(sum(Beta)+t(Beta)%*%(alpha))
   
   A=t(digamma(alpha.star)*(1+alpha)) ### 1xp
   B=digamma(rowSums(y)+alpha.star)%*%t(1+alpha) # nxp
   C=t(digamma(alpha.star.v)*(1+alpha)) ### 1xp
   D=digamma(y+matrix(alpha.star.v,numobs,p,byrow=TRUE))*matrix(1+alpha,numobs,p,byrow=TRUE) ### nxp
   
   out<- matrix(A,numobs,p,byrow=TRUE) - B - matrix(C,numobs,p,byrow=TRUE) + D
   return(out)} 
 
 
# alpha: funz. da massimizzare ####
alpha.f<-function(alpha,y,Beta,f.z.y){
  temp=calc.f.y.z(y,alpha,Beta)
  return(-sum(temp*f.z.y))
}  
 
 ## Funzioni per stima di Alpha #####
 
 # alpha: funz. per il gradiente #####
 alpha.f.grad<-function(alpha,y,Beta,f.z.y){
   p=ncol(y)
   numobs=nrow(y)
   temp1=matrix(f.z.y,numobs,p)*calc2grad.a(y,alpha,Beta)
   #return(list(grad=-apply(temp1,2,sum),hess=t(temp1)%*%(temp1)))
   return(-apply(temp1,2,sum))
 } 
 
 # alpha: calcolo gradiente #####
 calc2grad.a<-function(y,alpha,Beta){
   numobs=nrow(y)
   p=ncol(y)
   alpha.star.v=Beta*(1+alpha)
   alpha.star=c(sum(Beta)+t(Beta)%*%(alpha))
   
   A=t(digamma(alpha.star)*(Beta)) ### 1xp
   B=digamma(rowSums(y)+alpha.star)%*%t(Beta) # nxp
   C=t(digamma(alpha.star.v)*(Beta)) ### 1xp
   D=digamma(y+matrix(alpha.star.v,numobs,p,byrow=TRUE))*matrix(Beta,numobs,p,byrow=TRUE) ### nxp
   
   out<- matrix(A,numobs,p,byrow=TRUE) - B - matrix(C,numobs,p,byrow=TRUE) + D
   return(out)} 
 
 


###########################################################################

misc <-function(classification, truth){
    q <- function(map, len, x) {
      x <- as.character(x)
      map <- lapply(map, as.character)
      y <- sapply(map, function(x) x[1])
      best <- y != x
      if (all(len) == 1) 
        return(best)
      errmin <- sum(as.numeric(best))
      z <- sapply(map, function(x) x[length(x)])
      mask <- len != 1
      counter <- rep(0, length(len))
      k <- sum(as.numeric(mask))
      j <- 0
      while (y != z) {
        i <- k - j
        m <- mask[i]
        counter[m] <- (counter[m]%%len[m]) + 1
        y[x == names(map)[m]] <- map[[m]][counter[m]]
        temp <- y != x
        err <- sum(as.numeric(temp))
        if (err < errmin) {
          errmin <- err
          best <- temp
        }
        j <- (j + 1)%%k
      }
      best
    }
    if (any(isNA <- is.na(classification))) {
      classification <- as.character(classification)
      nachar <- paste(unique(classification[!isNA]), collapse = "")
      classification[isNA] <- nachar
    }
    MAP <- mapClass(classification, truth)
    len <- sapply(MAP[[1]], length)
    if (all(len) == 1) {
      CtoT <- unlist(MAP[[1]])
      I <- match(as.character(classification), names(CtoT), 
                 nomatch = 0)
      one <- CtoT[I] != truth
    }
    else {
      one <- q(MAP[[1]], len, truth)
    }
    len <- sapply(MAP[[2]], length)
    if (all(len) == 1) {
      TtoC <- unlist(MAP[[2]])
      I <- match(as.character(truth), names(TtoC), nomatch = 0)
      two <- TtoC[I] != classification
    }
    else {
      two <- q(MAP[[2]], len, classification)
    }
    err <- if (sum(as.numeric(one)) > sum(as.numeric(two))) 
      as.vector(one)
    else as.vector(two)
    bad <- seq(along = classification)[err]
    errorRate = length(bad)/length(truth)
    return(errorRate)
  }


is.odd=function(x) (x%%2==1)

norm_vec <- function(x) sqrt(sum(x^2))
norm_alpha<-function(x) 1.998*(x-min(x))/(max(x)-min(x))-0.999

