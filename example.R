## This is a code file of the first simulation (multimodal dataset simulation) in the RP-PCA paper.

method=c("CPCA", "T_PCA","S_ROB", "ROB","RP_PCA")
short_method=c("S_ROB", "ROB","RP_PCA")
measure=c("maxsub", "f.angle" ) #here we use these two measures to understand the performance of each method.


n.sim=20 #our simulation gives mean and standard deviations of 20 results.


main_a=0.75
jump=0.9

mult.m=list()
mult.jump=list()
for(dn in 1:4){mult.m[[dn]]=mult.jump[[dn]]=list()}
data_names=c("n-100-20","t-100-20","n-50-100","t-50-100")
names(mult.m)=names(mult.jump)=data_names

set.seed(1)
out.type.vector=c("Orthogonal", "Bad")[rbinom(n.sim,1,0.5)+1] #randomly deside that which data will be an orthogonal or a bad outlier set.
x.vector=runif(n.sim,10, 20) #this is the outlier center of each data set.



datalist=list()
# Four simulation data sets are generated for each seed number from 1 to 20.
# For example, the n-100-20 data of seed=3 is like the plots:
ex=simulation.data2(n=100, p=20, rank=5,
                 error.scale=0.05^2,
                 random.V = T,
                 contam.k=1, seed=3,
                 out.type="Bad",out.rate=0.15,  
                 out.parameter = x.vector[3] ,
                 score.sigma=diag(c(17,13.5,8,3,1)),data.plot = T )
#black lines denote data points before contamination (before error adding and outlier replacing), while red lines denote data points after contamination.

########################## simulation codes ##############################
for(seed in 1: n.sim){
  out.type=out.type.vector[seed]
  ## n-100-20 data generating
  datalist[[1]]=simulation.data2(n=100, p=20, rank=5,
                                 error.scale=0.05^2,
                                 random.V = T,
                                 contam.k=1, seed=seed,
                                 out.type=out.type,out.rate=0.15,  
                                 out.parameter = x.vector[seed] ,
                                 score.sigma=diag(c(17,13.5,8,3,1)))
  ## t-100-20 data generating
  datalist[[2]]=simulation.data2(n=100, p=20, rank=5,error.df=5 ,score.df = 5,
                                 error.scale=0.05^2,
                                 random.V = T,
                                 contam.k=1, seed=seed,
                                 out.type=out.type,out.rate=0.15,  
                                 out.parameter = x.vector[seed] ,
                                 score.sigma=diag(c(17,13.5,8,3,1)))
  #n-100-50 data generating
  datalist[[3]]=simulation.data2(n=50, p=100, rank=5,
                                 error.scale=0.05^2,
                                 random.V = T,
                                 contam.k=1, seed=seed,
                                 out.type=out.type,out.rate=0.15,  
                                 out.parameter = x.vector[seed] ,
                                 score.sigma=diag(c(17,13.5,8,3,1)))
  
  #t-100-50 data generating
  datalist[[4]]=simulation.data2(n=50, p=100, rank=5,error.df=5 ,score.df = 5,
                                 error.scale=0.05^2,
                                 random.V = T,
                                 contam.k=1, seed=seed,
                                 out.type=out.type,out.rate=0.15,  
                                 out.parameter = x.vector[seed] ,
                                 score.sigma=diag(c(17,13.5,8,3,1)))
  
  for(dn in 1:4){  
    
    # we use working manchine to simulate. 
    mult.m[[dn]][[seed]]= working.machine(data.object=datalist[[dn]],
                                          alpha=main_a, required.rank=5,
                                          method=method, measure=measure)$result
    temp_low=mult.m[[dn]][[seed]] ## this is the result of alpha=0.75 for j.f.angle or j.maxsub.
    temp_up= working.machine(data.object=datalist[[dn]],alpha=jump, required.rank=5,method=short_method, measure=measure)$result ## this is the result of alpha=0.9 for j.f.angle or j.maxsub.
    mult.jump[[dn]][[seed]]=list()
    for(msr in 1: length(measure)) mult.jump[[dn]][[seed]][[msr]]= temp_up[[msr]]-temp_low[[msr]][3:5] ##getting j.maxsub and j.f.angle.
    names(mult.jump[[dn]][[seed]])=paste("j.",measure, sep="")
    cat(paste("seed=", seed, "of set", data_names[dn],"\n" ))
   }
  
}



##this function sumarizes the results giving mean and standard deviation.
summary=function(types, method){
  joint=sds=means=list()
  for(msr in 1: length(measure)){
    O=t(sapply(types, function(Y) sapply(Y, function(X) X)[,msr]))
    colnames(O)=method
    means[[msr]]=round(apply(O, 2, mean),3)
    sds[[msr]]=round(apply(O, 2, sd),3)
    joint[[msr]]=paste(means[[msr]],"(", sds[[msr]], ")", sep = "")
  }
  names(joint)=names(sds)=names(means)=measure
  for(msr in 1: length(measure))  names(joint[[msr]])=method
  return(list(joint=t(sapply(joint, function(X)X)),means=t(sapply(means, function(X)X)),sds=t(sapply(sds, function(X)X))))
} 

##getting the summary results for four kinds of sets : "n-100-20","t-100-20","n-50-100", and "t-50-100"
mult.m.result=mult.jump.result=list();
for(dn in 1:4){
  mult.m.result[[dn]]=summary(mult.m[[dn]], method)
  mult.jump.result[[dn]]=summary(mult.jump[[dn]], short_method)
}

for(dn in 1:4){
cat(paste("\n ####################################################### \n\n","The results of dataset", data_names[dn],"\n","\n"))
print(mult.m.result[[dn]]$joint) 
print(mult.jump.result[[dn]]$joint)
cat("\n ####################################################### \n")
}
 
