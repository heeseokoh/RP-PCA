## Those are libraries for PCA and T-PCA
require(matrixcalc)
require(car)
require(rotations)
require(pracma)
## Those are libraries for ROBPCA and S-ROB
require(rrcov)
require(rospca)
## This is a library for ROBPCA 
require(robustbase)


###############################################################################
#######   CPCA and T-PCA and their related function   ######################### 
###############################################################################
 
#This is sign flip function called in CPCA(PCA) and T_PCA.
.signflip=function(V){
  return(apply(V, 2L, function(x) if (x[which.max(abs(x))] < 0) -x else x))
}

#This is function of conventional PCA.
CPCA<-PCA<-function(data,rank=dim(data)[2], signflip=T, ...){
  
  svd=svd(var(data), nv=rank, nu=rank)
  D=svd$d
  V=svd$v
  if(signflip) V=.signflip(V)
  mu=apply(data,2,mean)
  score=(crossprod(V,t(data)-mu))
  fit=V%*%score+mu
  obj=list(D=D,V=V,m=mu,fit=t(fit), score=t(score))
  return(obj)
}

### This is T-PCA function, and here nu is the key parameter!!! 
T_PCA<-function(Y, rank=(dim(Y)[2]-1),alpha=0.75, sigmasq0=0.01,nu0=c(2,0.5), screeplot=F, signflip=T, outlier.plot=F, n.direct=250, tol=0.00001, maxiter=1000, ridge=0, ...){
  
  #model : X=Wt+m+E
  
  ##input:
  #Y : nxd data, rank=q, explicit=F, tol, maxiter
 
  ##inner variables:
  #X: dxn data
  #i: 1, ..., d dimensions
  #j: 1, ..., n data numbers
  #m: MLE mean of data
  #Xbar: posterior mean of X |Y
  #Xsigma: posterior variance of X |Y
  #iter: total iteration number
  #N: total obervation number =n*d-length(omit)
  t.start=proc.time()
  iter=0
  X=t(Y)
  vacant.obs=which(apply(X,2, function(x) all(is.na(x))))
  
  if(!is.null(dim(vacant.obs))){
    cat(paste("\n will be disregarded vacant.obs =",vacant.obs ,sep=" " ))
    X=X[,-vacant.obs]
  }
  
  n=dim(X)[2]
  d=dim(X)[1]
  N=n*d
  
  #check completeness.
  omit=which(is.na((data.frame((X)))))
  N=n*d
  complete=(length(omit)==0)
  
  nu.diff<-W.diff<-rep(0,maxiter)
  
  ##the case of complete data
  if(complete){
    
    #initialize
    conduct=tryCatch(
      {
        pca0=ROB(Y, rank=rank, n.direct = n.direct,screeplot=screeplot, outlier.plot=outlier.plot)
        1
        
      }, error = function(e) {
        return(0)
      }, warning = function(w) {
        return(0)
      }
    )
    
    if (conduct==0){
      pca0=PCA(Y, rank=rank)
    }
    V=pca0$V
    m=pca0$m
    if(d>rank)  sigmasq=sum(diag(var(Y-pca0$fit)))/(d-rank) else sigmasq=sigmasq0
    
    for(trial in 1: length(nu0)){
      nu=nu0[trial]  
      ####end of initialize##
      
      ##iteration start## if there is an error, change nu, and try agian.
      conduct1=tryCatch(
        { 
          if(rank==1) W=V%*%(sqrt(apply(crossprod(V,(X-m)),1, var))) else W=V%*%diag(sqrt(apply(crossprod(V,(X-m)),1, var)))
          nu.past=nu
          W.past=W
          #two stage EM algorithm
          for(iter in 1: maxiter){
            
            #stage 1
            #Estep
            C=tcrossprod(W)+diag(rep(sigmasq, d))
            p=apply(hadamard.prod(t(crossprod(X-m, solve(C))),X-m),2, sum)
            u=(nu+d)/(p+nu)
            #Mstep
            m=as.vector(X%*%u/sum(u)) # m is updated at stage 1.
            
            #stage 2
            #Estep
            if(rank==1) M=crossprod(W)+ sigmasq  else  M=crossprod(W)+diag(rep(sigmasq, rank))
            p=apply(hadamard.prod(t(crossprod(X-m, solve(C))),X-m),2, sum)
            u=(nu+d)/(p+nu)
            M.inv=solve(M)
            t=tcrossprod(M.inv,W)%*%(X-m)
            ut=apply(t, 1, function(x) x*u)
            sum_utt=n*sigmasq*M.inv+t%*%ut
            logu=digamma((nu+d)/2)-log((nu+p)/2)
            nu=nlm(function(x) (-x/2*(mean(logu)+log(x/2)-mean(u))+log(gamma(x/2))), p=nu)$estimate
            #Mstep
            X.m=X-m
            W=(X.m)%*%ut%*%solve(sum_utt)
            sigmasq=as.numeric(( apply(hadamard.prod(X.m, X.m),2, sum)%*%u -2*sum(hadamard.prod(t(ut),crossprod(W,X.m)))+ sum(diag(crossprod(W)%*%sum_utt)))/N)
            
            #converge check
            W.diff[iter]=norm(W.past-W, "F")
            nu.diff[iter]=nu.past-nu
            
            if((W.diff[iter]<tol &nu.diff[iter]<tol) | iter>=maxiter){  
              break
            }else{
              W.past=W
            }#end of converge check and variable change to updated values
          }#end of iteration
          
          1},
        error = function(e) {
          return(0)
        })
      if(conduct1==1) break ##go to next step
      ##else, iterate agian with the next nu0
    }
  }else{
    ##EM algorithm of a missing data case.
    
    # understanding missing points
    d.vector=apply(X, 2, function(x)  sum(!is.na(x)))
    p=rep(0,n)
    iofj=list(); jofi=list(); imissofj=list()
    for(i in 1: d) {jofi[[i]]=which(! is.na(X[i,])); if(identical(jofi[[i]], integer(0))) stop(paste("remove only NA row:  the row of "), i)}
    for(j in 1: n) {iofj[[j]]=which(! is.na(X[,j]));
    imissofj[[j]]=which(is.na(X[,j]));
    if(identical(iofj[[j]], integer(0))) stop(paste("remove only NA row:  the col of "), j)
    }
    
    N=N-length(omit)
    
    #initialize
    set.seed(1)
    
    ## V, m: randomly generate
    V=rortho(d)[,1:rank]
    m=apply(Y, 2, mean, na.rm=T)
    
    sigmasq=sigmasq0
    
    for(trial in 1: length(nu0)){
      nu=nu0[trial]
      
      
      conduct1=tryCatch(
        {
          
          W=matrix(rnorm(d*rank), nrow=d, ncol=rank)
          sigmasq=runif(1, 0.1,1)
          
          ## end of the initialization.
          
          W.past=W
          ##iteration start## if there is an error, change nu, and try agian.
          #two stage EM algorithm
          for(iter in 1: maxiter){
            #stage 1
            #Estep
            C=(tcrossprod(W)+diag(d)*sigmasq)
            C.inv=solve(C)
            Z=X
            X.m=X-m
            for(j in 1:n){
              if(length(imissofj[[j]])==0){ #it means that jth observe does not have a missing
                p[j]=X.m[,j]%*%C.inv%*%X.m[,j]
              }else{
                Coo.inv=solve(C[iofj[[j]] ,iofj[[j]] ])
                xo.mo=X.m[iofj[[j]],j]
                coo.inv.xo.mo=Coo.inv%*%xo.mo
                p[j]=xo.mo%*%coo.inv.xo.mo
                Z[imissofj[[j]],j]=m[imissofj[[j]]]+C[imissofj[[j]],iofj[[j]]]%*%coo.inv.xo.mo
              }
            }
            u=(nu+d.vector)/(p+nu)
            #Mstep
            m=as.vector(Z%*%u/sum(u)) # m is updated at stage 1.
            
            #stage 2
            #Estep
            if(rank==1) M=crossprod(W)+ sigmasq  else  M=crossprod(W)+diag(rep(sigmasq, rank))
            X.m=X-m
            for(j in 1:n){
              if(length(imissofj[[j]])==0){ #it means that jth observe does not have a missing
                p[j]=X.m[,j]%*%C.inv%*%X.m[,j]
              }else{
                Coo.inv=solve(C[iofj[[j]] ,iofj[[j]] ])
                xo.mo=X.m[iofj[[j]],j]
                coo.inv.xo.mo=Coo.inv%*%xo.mo
                p[j]=xo.mo%*%coo.inv.xo.mo
                Z[imissofj[[j]],j]=m[imissofj[[j]]]+C[imissofj[[j]],iofj[[j]]]%*%coo.inv.xo.mo
              }
            }
            u=(nu+d.vector)/(p+nu)
            logu=digamma((nu+d.vector)/2)-log((nu+p)/2)
            Z.m=Z-m
            A=(Z.m)%*%apply(Z.m, 1, function(x) x*u)
            for(j in 1: n){
              if(length(imissofj[[j]])!=0){ #it means that jth observe has missings
                A[imissofj[[j]],imissofj[[j]]]=A[imissofj[[j]],imissofj[[j]]]+C[imissofj[[j]],imissofj[[j]]]-tcrossprod(C[imissofj[[j]],iofj[[j]]]%*%solve(C[iofj[[j]] ,iofj[[j]] ]),C[imissofj[[j]],iofj[[j]]])
              }
            }
            M.inv=solve(M)
            W.M.inv=W%*%M.inv
            sum_utt=n*sigmasq*M.inv+crossprod(W.M.inv,A)%*%W.M.inv
            sum_uxt=A%*%W.M.inv
            #Mstep
            W=sum_uxt%*%solve(sum_utt)
            sigmasq=(sum(diag(A))+sum(diag(crossprod(W)%*%sum_utt))-2*sum(diag(tcrossprod(sum_uxt, W)))) /(n*d)
            
            nu=nlm(function(x) (-x/2*(mean(logu)+log(x/2)-mean(u))+log(gamma(x/2))), p=nu)$estimate
            
            #convergence check
            W.diff[iter]=norm(W.past-W, "F")
            if(W.diff[iter]<tol | iter>=maxiter){  
              break
            }else{
              W.past=W
            }#end of convergence check and variable change to updated values
          }#end of iteration
          1
        },
        error = function(e) {
          return(0)
        })
      if(conduct1==1) break ##go to next step
      ##else, iterate agian with the next nu0
    }
    
    t=crossprod( W.M.inv,Z.m)
  }
  
  svd=svd(tcrossprod(W), rank, rank)
  V=svd$v
  if(signflip) V=.signflip(V)
  
  Y.hat=W%*%solve(crossprod(W))%*%M%*%t+m
  if(!is.null(dim(vacant.obs))){
    Y[-vacant.obs, ]=t(Y.hat)
  }else{
    Y=t(Y.hat)
  }
  score=crossprod(V,Y.hat-m)
  t.finish=proc.time()
  
  D=svd(var(t(Y.hat)))$d[1:rank]
  
  return(list(fit=Y, m=m , V=V,D=D, Yvar=tcrossprod(W)+diag(d)*sigmasq,  sigmasq=sigmasq , iter=iter,time=t.finish-t.start, w.diff=W.diff[1:iter], score=t(score), nu=nu))
  
}




###############################################################################
#######   functions for ROBPCA and R-ROB    ################################## 
###############################################################################

# ROB is function just get information PcaHubert implemented in rrcov, and through ROB the information has the form of what is suitable
# to our working machine codes. 
# Here, mcd=F because otherwise, PcaHubert can output not ROBPCA results but just MCD results.
ROB<-function(Y, alpha=0.75, rank=0, screeplot=T, outlier.plot=T, robust.scaling=F, inform.rank=F,n.direct=250,...){
  #This function will not change n.direct. Therefore, the argument n.direct is a fake variable.
  rob=PcaHubert(Y, k=rank, mcd=F, alpha=alpha, scale=robust.scaling,signflip = T)
  if(screeplot) screeplot(rob)
  if(outlier.plot) {
    plot(rob, main=paste("ROB, No. of outliers= ",length(which(!rob@flag) ) ));
  }
  V=rob@loadings
  m=rob@center
  D=rob@eigenvalues
  
  # to get fitted values which are not provided by PcaHubert
  tYc=t(Y)-m
  score.c=tcrossprod(V)%*%tYc
  fit= t(score.c + m)
  return(list(D=D,V=V,m=m,fit=fit, rob.obj=rob,score=rob@scores, flag=which(!rob@flag)))
}

# S_ROB is function just get information robpca with skew=1 (that is, S-ROB result) implemented in rrcov, 
# and through S_ROB the information has the form of what is suitable to our working machine codes. 
# Here, mcd=F because otherwise, robpca might output not S_ROB results but just MCD results.

S_ROB<-function(Y, alpha=0.75, rank=0, screeplot=T, outlier.plot=T, robust.scaling=F, inform.rank=F,n.direct=250,...){
  #This function will not change n.direct.
  
  rob=robpca(Y, k=rank, mcd=F,skew=T, alpha=alpha, scale=robust.scaling,signflip = T)
  
  V=rob$loadings
  m=rob$center
  D=rob$eigenvalues
  
  # to get fitted values which are not provided by robpca
  tYc=t(Y)-m
  score.c=tcrossprod(V)%*%tYc
  fit= t(score.c + m)
  
  return(list(D=D,V=V,m=m,fit=fit, rob.obj=rob,score=rob$scores, flag=which(!rob$flag.all)))
}


###############################################################################
#######   RP-PCA and itsrelated functions   ################################## 
###############################################################################

# RP-PCA function exploits PCA and T_PCA as well as g functions and others below.

# This is sign flip function used for ordering PC directions in RP-PCA.
# input :: V: PC vectors, t: scores, rank: rank.
# output :: ordinated V and t.
   .multi_signflip=function(V, t, rank=dim(V)[2]){
    for(i in 1: rank){
      x=V[,i]
      if(x[which.max(abs(x))] < 0){
        V[,i]=-x
        t[i,]=-t[i,]
      }
    }
    return(list(V=V, t=t))
  }
  
  # This is function giving outl of each data observation in a data set.
  OUTL= function(t, alpha=0.75){
    x=covMcd(t, alpha=alpha)
    if(x$cov>10^(-12)) return( abs(t-x$center)/x$cov)
    else return(t*0)
  }
  
  # Gp, GBOX and God are g functions which takes data and gives resultant g-vector composed of 0 or 1. 
  Gp=function(data, alpha=0.75, n.direct=200){
    X=t(data)
    p=dim(X)[1]
    n=dim(X)[2]
    h=ceiling(alpha*n)
    set.seed(251120134)
    
    
    random=matrix(ceiling(runif(2*n.direct, min=1, max=n)), nrow=2)
    B=X[,random[1,]]-X[,random[2,]] # p x 250 (or p x n.direct)
    Bnorm=apply(B, 2, norm, "2")
    B=B[,which(Bnorm>10^(-12))]
    Bnorm=Bnorm[which(Bnorm>10^(-12))]
    m=length(Bnorm)
    A=(apply(B, 1, function(x) x/Bnorm) ) #length 1 ramdom(?) direction vectors.
    Y=A%*%X
    
    Z=apply(Y, 1, OUTL, alpha=alpha)
    d=apply(Z, 1, max)
    
    g_vector=(d<=sort(d)[h])
    
    return(g_vector)
  }
  GBOX=function(data, alpha=0.75, n.direct=200){
    
    X=t(data)
    p=dim(X)[1]
    n=dim(X)[2]
    h=ceiling(alpha*n)
    set.seed(251120134)
    
    random=matrix(ceiling(runif(2*n.direct, min=1, max=n)), nrow=2)
    B=X[,random[1,]]-X[,random[2,]] # p x 250 (or p x n.direct)
    Bnorm=apply(B, 2, norm, "2")
    B=B[,which(Bnorm>10^(-12))]
    Bnorm=Bnorm[which(Bnorm>10^(-12))]
    m=length(Bnorm)
    A=(apply(B, 1, function(x) x/Bnorm) ) #length 1 ramdom(?) direction vectors.
    Y=A%*%X
    
    Z=apply(Y, 1, OUTL, alpha=alpha)
    d=apply(Z, 1, max)
    
    gp_vector=(d<=sort(d)[h])
    
    su=base::summary(d[gp_vector])
    whisker=2.5*su[5]-1.5*su[2]  ##This is q3+1.5IQR (upper whisker of boxplot of h points)
    
    g_vector=( gp_vector | (d<whisker)) 
    return(g_vector)
  }
  God=function(data, rank=dim(data)[2]-1,alpha=0.75, n.direct=250){
    unit=t_robust_unit("Gp",  rank=rank, g_options=list(data=data,alpha=alpha,n.direct=n.direct) )
    
    m=unit$pca$m
    V=unit$pca$V
    p=length(m)
    Y=t(data)
    
    OD=apply(((diag(p)-tcrossprod(V))%*%(Y-m)), 2, norm, "2")
    mcd=covMcd(OD^(2/3), alpha=alpha)
    mu=mcd$center
    sigma=sqrt(mcd$cov)[1,1]
    COD=(mu+sigma*qnorm(0.975))^(3/2)
    g_vector=(OD<=COD)
    return(g_vector)
  }
  
  # t_robust_unit is function of t-robust unit in RP-PCA paper.
  # This is the key basic operator of RP-PCA.
  t_robust_unit=function(g="Gp",rank,g_options,tol= 1e-05, maxiter=1000){
    dim=dim(g_options$data)
    missing=any(is.na(g_options$data))
    p=dim[2]
    n=dim[1]
    iter=NA
    
    if(!missing){
      g_vector=do.call(g, c(g_options))
      
      pca=T_PCA(g_options$data[g_vector,], rank=rank,tol=tol, maxiter=maxiter)
      
      outs=g_options$data[!g_vector,]
    }else{
      keep=g_options$data
      Y=keep ##Y would be an imputed dataset.
      
      #1.prelimarily imputation via any mathod
      missloc=which(is.na(Y))
      pp=PPCA(Y, rank=rank)
      Y[missloc]=pp$fit[missloc]
      
      g_options$data=Y #impute
      
      g_vector=do.call(g, c(g_options)) #getting filter for primarily imputed data
      
      for(iter in 1:maxiter){
        g_vector_paste=g_vector
        
        g_missloc=which(is.na(keep[g_vector,]))
        
        not_g_missloc=which(is.na(keep[!g_vector,]))
        
        pca=T_PCA(Y[g_vector,], rank=rank)
        
        ###impute process###
        Y[g_vector,][g_missloc]=pca$fit[g_missloc] #impute
        ###### treat filtered outliers (also impute based on pca results) #### 
        if(sum(!g_vector)>0){
          outs=keep[!g_vector,]
          C=pca$Yvar
          m=pca$m
          if(is.null(dim(outs))) {
            if (any(is.na(outs))){
              nacomp=is.na(outs)
              outs[nacomp]=C[nacomp,!nacomp]%*%solve(C[!nacomp,!nacomp])%*%(outs[!nacomp]-m[!nacomp])
              
            }
          }else{
            for(i in which(apply(outs, 1, function(x) any(is.na(x))))){
              nacomp=is.na(outs[i,])
              outs[i,][nacomp]=C[nacomp,!nacomp]%*%solve(C[!nacomp,!nacomp])%*%(outs[i,][!nacomp]-m[!nacomp])
            }
          }
          Y[!g_vector,]=outs ###### impute for outliers
        }
        ##get a new filter
        g_options$data=Y
        g_vector=do.call(g, c(g_options))
        if(all(g_vector==g_vector_paste)) break
      }
    }
    ##fit and scores for outliers
    
    if(sum(!g_vector)>0){
      fit=matrix(0, ncol=p, nrow=n)
      score=matrix(0, ncol=rank, nrow=n)
      V=pca$V
      m=pca$m
      score[g_vector,]=pca$score
      
      if(is.null(dim(outs))){
        outscore=crossprod(V,((outs)-m))
      }else{
        outscore=crossprod(V,(t(outs)-m))
      }
      
      score[!g_vector,]=t(outscore)
      
      fit[g_vector,]=pca$fit
      fit[!g_vector,]= t(V%*%outscore+m)
      
      pca$fit=fit
      pca$score=score}
    return(list(pca=pca, g_vector=g_vector, iter=iter)) #iteration default is NA : no missing!
  }
  
  
  # This is RP-PCA function
  RP_PCA=function(data, rank=dim(data)[2]-1,alpha=0.75, n.direct=250,
                  pickup.rate=0, signflip=T, tol= 1e-05, maxiter=1000,...){
    
    W=t(data)
    p=dim(W)[1]
    n=dim(W)[2]
    
    complete=!any(is.na(data))
    
    ##rotate##
    #### Stage 0 ->rotation  
    if(complete){
      svd0=svd(var(t(W)))
      k0=max(which(svd0$d>10^(-12)))
      V_W=svd0$v[, 1:k0]
      mu_W=apply(W, 1, mean)
      Y=crossprod(V_W, W-mu_W)  # : k0 x n
      
      center=mu_W # : p x 1
      V=V_W # : p x k0
    }else{
      Y=W
      center<-mu_W<-rep(0,p)
      V<-V_W<-diag(p)#
    }
    
    ###stage 1;
    
    pca1=t_robust_unit(g="God",rank=rank,g_options=list(data=t(Y),rank=rank, alpha=alpha,n.direct=n.direct), tol=tol, maxiter=maxiter)
    mu_Y=pca1$pca$m
    V_Y=pca1$pca$V
    
    ###stage 2;
    if(rank>1){
      pca2=t_robust_unit(g="GBOX", rank=rank, g_options=list(data=pca1$pca$score,alpha=alpha,n.direct=n.direct), tol=tol , maxiter=maxiter)
      mu_X=pca2$pca$m
      V_X=pca2$pca$V
      
    }else{ 
      #if rank=1; there is no way to improve PC vectors 
      #the thing we can do is only re-centering
      #But, here, we do not anything. 
      #if someone want to improve this also, one might use univariate mcd like{mcd=covMcd(pca1$pca$score, alpha = alpha); mu_X=as.vector(mcd$center)}
      pca2=NULL;
      mu_X=0
      V_X=1
    }
    
    center=as.vector(center+V%*%mu_Y) # p x 1
    V=V%*%V_Y # p x rank
    
    center=as.vector(center+V%*%mu_X)
    V=V%*%V_X
    
    
    if(is.na(pca1$iter)){
      score=t(crossprod(V,W-center)) ## Do this because scores of t-pca results tend to be shrinked
    }else{
      score=pca1$pca$score  ## Do this because there are missing in W.
    }
    
    
    fit=t(V_W%*%t(pca1$pca$fit)+mu_W)
    D=pca1$pca$D
    
    if(signflip) {
      flip= .multi_signflip(V, score,rank)
      V=flip$V
      score=flip$t
    }
    
    
    return(list(fit=(fit),iter1=pca1$iter,iter2=pca2$iter,score=(score), m=center, V=V, D=D,flag1=pca1$g_vector, flag2=pca2$g_vector, nu1=pca1$pca$nu, nu2=pca2$pca$nu, u1=pca1$pca$u, u2=pca2$pca$u))
  } 
  # iter2 of RP-PCA result is always NA when there are missing values!! Because the Stage 2 does not impute to any points!
  # if there are no missing value, both iter1 and iter 2 are NA. 
  # iter1 and iter2 just the iteration number of imputation. 
 
