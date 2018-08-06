library(coda)

#### auxiliary functions common to both mcmc sampling ####

# values of lambda_n and a_n (see appendix 1)
# epsilon and beta are the ones described in the paper, N is the tip number
lambda_N=function(epsilon,beta,N,upper=10000){c(0,sapply(2:N, function(x){lambda.epsilon(epsilon,beta,x,upper)}))}
a_N=function(beta,N){c(0,sapply(2:N,function(n){sum(sapply(1:(n-1),function(j){beta(beta+1+j,beta+1+n-j)/((n+1)*beta(1+j,1+n-j))}))}))}

# probability of the sampled node position knowing the tree topology, the number of 
# unsampled nodes and the interval sizes
aux_lik=function(M,beta,alpha){
  n=nrow(M)
  int=M[,"int"]
  depths=M[,"depth"]
  proba=0
  current=c(1)
  
  while (length(current)>0){
    i=which.max(depths[current])
    p=alpha*(int[current[i]])-log(sum(exp(int[current]*alpha)))
    proba=p+proba
    
    
    current=insert(current,
                   (1:n)[M[,"father"]==M[current[i],"node"]],
                   i,T)
  }
  return(proba)
}

# compute the splits at each node of a tree
split=function(tree){
  n=tree$Nnode+1
  rep=matrix(0,3,2*n-1)
  rep[1,1:n]=rep(1,n)
  for (i in (2*n-2):1){
    rep[1,tree$edge[i,1]]=rep[1,tree$edge[i,1]]+rep[1,tree$edge[i,2]]
    if (rep[3,tree$edge[i,1]]==0){
      rep[3,tree$edge[i,1]]=rep[1,tree$edge[i,2]]
    }else{
      rep[2,tree$edge[i,1]]=rep[1,tree$edge[i,2]]
    }
  }
  return(rep)
}

insert=function(v,e,i,replace){
  if (i==0){
    c(e,v)
  }else  if(i==1){
    if (replace){
      c(e,v[-1])
    }else{c(v[1],e,v[-1])}
  }else if (i==length(v)){
    if (replace){
      c(v[-length(v)],e)
    }else{
      c(v,e)
    }
  }else if (replace){
    c(v[1:(i-1)],e,v[(i+1):length(v)])
  }else{
    c(v[1:i],e,v[(i+1):length(v)])
  }
}

rbactrian <- function(n, m=0.95){
  mBactrian = m
  sBactrian = sqrt(1-m^2)
  z = mBactrian + rnorm(n,0,1)*sBactrian
  rdunif <- runif(n,0,1)<0.5
  sign <- ifelse(rdunif,-1,1)
  z=z*sign
  return(z)
}

#### inference function of alpha only ####

# 3 functions that are used within the inference function :

# write the tree in a table
transform=function(tree,unsampled,list_depth=NULL,change_depth=0){
  split=split(tree)
  ntip=tree$Nnode+1
  depth=c(rep(0,ntip),depth(tree))
  p_unsampled=0
  compute_depth=is.null(list_depth)
  if(compute_depth){list_depth=list(c())}
  fils=tree$edge[tree$edge[,1]==(ntip+1),2]
  M=c(ntip+1,-1,depth[ntip+1],ntip)
  k=2*ntip
  id.unsampled=1
  for (i in (ntip+2):(2*ntip-1)){
    father=tree$edge[tree$edge[,2]==i,1]
    node_father=father
    pf=depth[i]
    pp=depth[father]
    fils=tree$edge[tree$edge[,1]==i,2]
    if (unsampled[i]>0){
      if(compute_depth | (change_depth==i)){
        depths=runif(unsampled[i],pf,pp)
        depths=depths[order(depths,decreasing = T)]
        list_depth[[i]]=depths
      }else{
        depths=list_depth[[i]]
      }
      for (j in 1:unsampled[i]){
        M=rbind(M,c(k,father,depths[j],split[1,i]))
        father=k
        k=k+1
        id.unsampled=id.unsampled+1
      }
      p_unsampled=p_unsampled+sum(log(1:unsampled[i]))-log(pp-pf)*unsampled[i]
    }
    M=rbind(M,c(i,father,pf,split[1,i]))
  }
  colnames(M)=c("node","father","depth","ntip")
  for(j in 1:M[1,"depth"]){
    N=sum(M[,"depth"]>(j-1) & M[,"depth"]<(j))
    if(N>0){
      p_unsampled=p_unsampled-sum(log(1:N))
    }
  }
  
  return(list(M=M,p_unsampled=p_unsampled,list_depth=list_depth))
}

# add unsampled nodes and interval knowledge
enhance=function(tree,alpha,beta,epsilon,lambdaN,aN){
  N=tree$Nnode
  unsampled=c(-1,rep(-1,2*N))
  p_nUnsampled=0
  
  for (i in (N+3):(2*N+1)){
    n=extract.clade(tree,i)$Nnode+1
    
    if(n==1){unsampled[i]=-1
    }else {
      lambda=lambdaN[n]
      an=aN[n]
      sigma=rexp(1,an)
      unsampled[i]=rpois(1,lambda*sigma)
      p_nUnsampled=p_nUnsampled+log(an)-(unsampled[i]+1)*log(lambda+an)+(unsampled[i])*log(lambda)
      if(is.na(unsampled[i])) print(paste("a",alpha,"b",beta,"t",unsampled[i],"lamb",lambda,"an",an,"sig",sigma,"n",n))
    }
  }
  
  tr=transform(tree,unsampled)
  p_unsampled=tr$p_unsampled
  list_depth=tr$list_depth
  tr=tr$M
  
  int=c(0,rep(-Inf,nrow(tr)-1))
  Rs=c(1,rep(-1,nrow(tr)-1))
  
  for(i in 2:nrow(tr)){
    if(int[i]==-Inf){
      father=tr[i,2]
      j=(1:nrow(tr))[tr[,1]==father]
      if(father<(2*N+2)){
        R=simulate.R(beta,tr[j,4],tr[i,4])[1]
        int[i]=int[j]+log(R)
        Rs[i]=R
        if(i<nrow(tr)){
          k=((i+1):nrow(tr))[tr[(i+1):nrow(tr),2]==father]
          int[k]=int[j]+log(1-R)
          Rs[k]=(1-R)
        }
      }else{
        R=exp(-1*simulate.Yi(1,epsilon,beta,tr[i,4]))
        int[i]=int[j]+log(R)
        Rs[i]=R
      }
    }
  }
  
  
  tr=cbind(tr,int)
  return(list(transform=tr,unsampled=unsampled,int=int,Rs=Rs,p_nUnsampled=p_nUnsampled,p_unsampled=p_unsampled,list_depth=list_depth,Infinite=F))
}

# change the number of unsampled nodes and the intervals for a given node set
change_int=function(enhanced,ind,tree,alpha,beta,epsilon,lambdaN,aN){
  changeMax=5
  N=tree$Nnode
  unsampled=enhanced$unsampled
  int=enhanced$int
  p_nUnsampled=enhanced$p_nUnsampled
  list_depth=enhanced$list_depth
  continue=T
  if(is.null(ind) ){
    rm=c()
  }else if (ind==(N+2)) {
    rm=c()
  }else{
    n=extract.clade(tree,ind)$Nnode+1
    rm=c()
    tr=enhanced$transform
    node=which(tr[,1]==ind)
    if(unsampled[ind]>0){
      rm=(node-1):(node-unsampled[ind])
      node=node-unsampled[ind]
    }
    
    if(n==1){unsampled[ind]=-1
    }else if(ind==(N+2)){unsampled[ind]=-1}else{
      lambda=lambdaN[n]
      an=aN[n]
      p_nUnsampled=p_nUnsampled+(unsampled[ind]+1)*log(lambda+an)-(unsampled[ind])*log(lambda)
      u=runif(1)
      if(u<(0.4)){
        unsampled[ind]=unsampled[ind]-1
      }else if(u>0.6){
        unsampled[ind]=unsampled[ind]+1
      }
      unsampled[ind]=unsampled[ind]+floor(rnorm(1,mean = 0.5,sd = 10))
      if(unsampled[ind]<0){
        proba.unsampled=-Inf
        continue=F
      }else{
        p_nUnsampled=p_nUnsampled-(unsampled[ind]+1)*log(lambda+an)+(unsampled[ind])*log(lambda)
        continue=T
      }
      if(is.na(unsampled[ind])) print(paste("a",alpha,"b",beta,"t",unsampled[ind],"lamb",lambda,"an",an,"sig",sigma,"n",n))
    }
  }
  
  if(! continue){
    enhanced$p_nUnsampled=-Inf
    return(enhanced)
  }
  
  tr=transform(tree,unsampled,list_depth = list_depth,change_depth = max(ind,0))
  p_unsampled=tr$p_unsampled
  list_depth=tr$list_depth
  tr=tr$M
  oldRs=c()
  
  if(length(rm)>0){
    oldRs=enhanced$Rs[rm]
    Rs=enhanced$Rs[-rm]
    depth=enhanced$transform[-rm,"depth"]
    
  }else{
    
    Rs=enhanced$Rs
    depth=enhanced$transform[,"depth"]
    
  }
  if(is.null(ind)){
    Rs=enhanced$Rs
  }else{if(unsampled[ind]>-1){
    add=c(-1,sample(c(oldRs,rep(-1,changeMax+max(0,unsampled[ind]-length(oldRs)))),unsampled[ind]))
    if(node<length(Rs)){
      Rs=c(Rs[1:(node-1)],add,Rs[(node+1):length(Rs)])
      depth=c(depth[1:(node-1)],tr[node:(node+unsampled[ind]),"depth"],depth[(node+1):length(depth)])
      
    }else{
      Rs=c(Rs[1:(node-1)],add)
      depth=c(depth[1:(node-1)],tr[node:(node+unsampled[ind]),"depth"])}
  }}
  
  tr[,"depth"]=depth
  
  int=c(0,rep(-Inf,nrow(tr)-1))
  if(! is.null(ind)){if(ind>(N+2)){father=tr[node,2]
  Ks=which(tr[,2]==father)
  j=(1:nrow(tr))[tr[,1]==father]
  if(Rs[Ks[1]]==1){R=simulate.R(beta,tr[j,4],tr[Ks[1],4])[1]
  Rs[Ks[1]]=R
  if(length(Ks)==2){Rs[Ks[2]]=(1-R)}}}}
  Rs[1]=1
  
  for(i in 2:nrow(tr)){
    if(int[i]==-Inf){
      father=tr[i,2]
      j=(1:nrow(tr))[tr[,1]==father]
      if(father<(2*N+2)){
        if(Rs[i]==-1){
          R=simulate.R(beta,tr[j,4],tr[i,4])[1]
          Rs[i]=R}
        R=Rs[i]
        int[i]=int[j]+log(R)
        if(i<nrow(tr)){
          k=((i+1):nrow(tr))[tr[(i+1):nrow(tr),2]==father]
          int[k]=int[j]+log(1-R)
          Rs[k]=(1-R)
        }
      }else{
        if(Rs[i]==-1){
          R=exp(-1*simulate.Yi(1,epsilon,beta,tr[i,4]))
          Rs[i]=R}
        R=Rs[i]
        int[i]=int[j]+log(R)
        
      }
    }
  }
  if(any(Rs<0) ) {print(ind); debug(change_int); change_int(enhance,ind,tree,alpha,beta,epsilon)}
  tr=cbind(tr,int)
  return(list(transform=tr,unsampled=unsampled,int=int,Rs=Rs,p_nUnsampled=p_nUnsampled,p_unsampled=p_unsampled,list_depth=list_depth,Infinite=F))
}

# and the inference function :

# run the inference of alpha
# tree is the phylogeny, epsilon is the minimum split size (as described i appendix 1),
# beta is the MLE of the beta statistic for the tree, 
# lambadN and aN are the vector given by the functions lambda_N and a_N, 
# niter is the number of iterations in the mcmc chain,  
# ini is the initial value of alpha, V the initial scaling factor fot the proposal function,
# [ma,Ma] is the interval within which alpha is allowed to be, 
# verbose is the number of iteration after which the state of the mcmc is printed (if silent == F)
# the proposal is scaled every Nadapt iterations (if at least NadaptMin iterations have been performed)
# chain is either NULL (and a new inference is being performed) or the output of mcmc_alpha (and the chain is continued)
# An example of the use of this function is given in the end of this document

mcmc_alpha=function(tree,epsilon,beta,lambdaN,aN,niter,ini=0,V=0.1,chain=NULL,
                    verbose=10,silent=T,Nadapt=Inf,NadaptMin=10,
                    ma=-4,Ma=4,proposal="bactrian",accOpt=0.3,Vmin=0.001){
  
  
  ND=rank(nodes.depths.ordonnes(tree))
  tree2=build.tree(ND)
  tree$edge.length=tree2$edge.length
  N=tree$Nnode
  
  ntip=tree$Nnode+1
  ordre=rank(max(depth(tree))+1-depth(tree))
  X=rep(1,ntip-1)
  for (i in 1:(ntip-1)){X[ordre[i]]=((ntip+1):(ntip*2-1))[i]}
  
  
  posterior2=function(alpha,enhance){
    if(enhance$Infinite){return(-Inf)}
    Y=try(aux_lik(enhance$transform,beta,alpha))
    if(inherits(Y,"try-error")){Y=-Inf}
    p2=enhance$p_nUnsampled
    p3=enhance$p_unsampled
    if(any(is.na(c(Y,p2,p3))) | any(is.nan(c(Y,p2,p3))) | any(is.infinite(c(Y,p2,p3)))){return(-Inf)}
    return(Y-p3+p2)
  }
  
  posterior=function(alpha,enhance){
    if(enhance$Infinite ){return(-Inf)}
    Y=try(aux_lik(enhance$transform,beta,alpha))
    if(inherits(Y,"try-error")){Y=-Inf}
    if(is.na(Y) | is.nan(Y) | is.infinite(Y)){return(-Inf)}
    return(Y)
  }
  
  if (is.null(chain)){
    a=ini
    j=1
    post_actuel=-Inf
    while(post_actuel==-Inf){
      int=enhance(tree,a,beta,epsilon,lambdaN,aN)
      post_actuel=posterior(a,int)
      post_actuel_int=posterior2(a,int)
    }
    Mcmc=c(a,sum(int$unsampled[-(1:(N+2))]),post_actuel)
  }else{
    a=chain$a
    int=chain$int
    Mcmc=chain$mcmc
    j=nrow(chain$mcmc)
    V=chain$V
    post_actuel=chain$post_actuel
    post_actuel_int=posterior2(a,int)
  }
  n=0
  for(i in 1:niter){
    if((i+j) %% Nadapt ==0 & (i+j)>NadaptMin){ 
      if(proposal=="uniform" | proposal=="bactrian"){
        acc=sum(diff(Mcmc[-c(1:floor(nrow(Mcmc)-Nadapt)),1])>0)/(Nadapt-1)
        if(acc<0.001){
          V=V/100
        }else if(acc>0.999){
          V=V*100
        }else{
          V = V * (tan((pi/2)*acc)/tan((pi/2)*accOpt))
        }
        V=max(V,Vmin)
      }else{
        V=max(Vmin,sqrt(var(Mcmc[-c(1:floor(nrow(Mcmc)/2)),1]))*2.36)
      }
      if(is.na(V)) V=1
    }
    
    #On tire les paramÃ¨tres proposÃ©s : 
    post_actuel=try(posterior(a,int))
    if(proposal=="uniform"){
      new_a = a + (runif(1,0,1)-0.5)*3.4641016*V
    }else if (proposal=="bactrian"){
      new_a =a + rbactrian(1)*V
    }else{
      new_a=rnorm(1,a,V)}
    if(new_a>Ma| new_a<ma){
      post=-Inf
    }else{
      post=try(posterior(new_a,int))}
    
    if(post>=post_actuel){
      a=new_a
      post_actuel=post
      post_actuel_int=posterior2(a,int)
    }else{
      u=runif(1,0,1)
      if (log(u)<(post-post_actuel)){
        a=new_a
        post_actuel=post
        post_actuel_int=posterior2(a,int)
      }
    }
    
    ks=(N+3):(2*N+1)
    for(k in ks){
      new_int=change_int(int,k,tree,a,beta,epsilon,lambdaN,aN)
      post=posterior2(a,new_int)
      
      if(post>=post_actuel_int){
        int=new_int
        post_actuel_int=post
        post_actuel=posterior(a,int)
      }else{
        u=runif(1,0,1)
        if (log(u)<(post-post_actuel_int)){
          int=new_int
          post_actuel_int=post
          post_actuel=posterior(a,int)
        }
      }}
    
    
    Mcmc=rbind(Mcmc,c(a,sum(int$unsampled[-(1:(N+2))]),post_actuel))
    
    n=n+1
    if(n==verbose){
      n=0
      print(paste(i,"alpha",Mcmc[i+j,1],"beta",beta,"post",Mcmc[i+j,3]))
    }
  }
  colnames(Mcmc)=c("alpha","nUnsampled","log_post")
  return(list(mcmc=mcmc(Mcmc),post_actuel=post_actuel,a=a,int=int,tree=tree,V=V))
}


#### inference fuction of alpha and eta ####

# 3 functions used within the inference function :

# write the tree in a table
Rmin=log(1e-10)-1
maxChange=15
transform_eta=function(tree,unsampled,eta,list_depth=NULL,change_depth=0){
  split=split(tree)
  ntip=tree$Nnode+1
  depth=c(rep(0,ntip),depth(tree))
  A1=-1
  A2=-1
  p_unsampled=0
  compute_depth=is.null(list_depth)
  if(compute_depth){list_depth=list(c())}
  fils=tree$edge[tree$edge[,1]==(ntip+1),2]
  if(fils[1]<=ntip) A1=tree$tip.ab[fils[1]]
  if(fils[2]<=ntip) {
    if (A1==-1) {A1=tree$tip.ab[fils[2]]}else{A2=tree$tip.ab[fils[2]]}}
  M=c(ntip+1,-1,depth[ntip+1],ntip,A1,A2)
  k=2*ntip
  id.unsampled=1
  for (i in (ntip+2):(2*ntip-1)){
    A1=-1
    A2=-1
    father=tree$edge[tree$edge[,2]==i,1]
    node_father=father
    pf=depth[i]
    pp=depth[father]
    fils=tree$edge[tree$edge[,1]==i,2]
    if(fils[1]<=ntip) A1=tree$tip.ab[fils[1]]
    if(fils[2]<=ntip) {
      if (A1==-1) {A1=tree$tip.ab[fils[2]]}else{A2=tree$tip.ab[fils[2]]}}
    if (unsampled[i]>0){
      if(compute_depth | (change_depth==i)){
        depths=runif(unsampled[i],pf,pp)
        depths=depths[order(depths,decreasing = T)]
        list_depth[[i]]=depths
      }else{
        depths=list_depth[[i]]
      }
      for (j in 1:unsampled[i]){
        M=rbind(M,c(k,father,depths[j],split[1,i],-1,-1))
        father=k
        k=k+1
        id.unsampled=id.unsampled+1
      }
      p_unsampled=p_unsampled+sum(log(1:unsampled[i]))-log(pp-pf)*unsampled[i]
    }
    M=rbind(M,c(i,father,pf,split[1,i],A1,A2))
  }
  
  colnames(M)=c("node","father","depth","ntip","A1","A2")
  for(j in 1:M[1,"depth"]){
    N=sum(M[,"depth"]>(j-1) & M[,"depth"]<(j))
    if(N>0){
      p_unsampled=p_unsampled-sum(log(1:N))
    }
  }
  return(list(M=M,p_unsampled=p_unsampled,list_depth=list_depth))
}

# add unsampled nodes and interval knowledge
enhance_eta=function(tree,alpha,beta,eta,epsilon,lambdaN,aN,uns){
  N=tree$Nnode
  unsampled=c(-1,rep(-1,2*N))
  p_nUnsampled=0
  
  for (i in (N+2):(2*N+1)){
    n=extract.clade(tree,i)$Nnode+1
    
    if(n==1){unsampled[i]=-1
    }else {
      lambda=lambdaN[n]
      an=aN[n]
      sigma=rexp(1,an)
      unsampled[i]=rpois(1,lambda*sigma)
      p_nUnsampled=p_nUnsampled+log(an)-(unsampled[i]+1)*log(lambda+an)+(unsampled[i])*log(lambda)
      if(is.na(unsampled[i])) print(paste("a",alpha,"b",beta,"t",unsampled[i],"lamb",lambda,"an",an,"sig",sigma,"n",n))
    }
  }
  
  tr=transform_eta(tree,unsampled,eta)
  p_unsampled=tr$p_unsampled
  list_depth=tr$list_depth
  tr=tr$M
  
  int=c(-0,rep(-Inf,nrow(tr)-1))
  Rs=c(0,rep(-Inf,nrow(tr)-1))
  As=rep(-1,nrow(tr))
  
  for(i in 1:N){
    do=(tr[,2]>=(2*N+2) & Rs==-Inf & tr[,4]==(i+1))
    if(sum(do)>0){Rs[do]=uns[[i]][sample(10000,sum(do),replace=T)]}
  }
  
  nTr= nrow(tr)
  nodes=rep(T,nTr)
  depth=tr[,"depth"]
  while(sum(nodes)>0){
    i=which.min(depth)
    node=tr[i,"node"]
    
    if(!((tr[i,"A1"]==-1) | (tr[i,"A2"]==-1 & tr[i,"node"]<(2*N+2)))){
      depth[i]=Inf
      if(tr[i,"A2"]>-1){
        As[i]=tr[i,"A1"]+tr[i,"A2"]
      }else{
        k=which(tr[,"father"]==node)
        As[i]=tr[i,"A1"]*(1+exp(eta*(log(1-exp(Rs[k]))-(Rs[k]))))
      }
      father=tr[i,"father"]
      nodes[i]=F
      if(father>-1) {
        j=which(tr[,"node"]==father)
        if(tr[j,"A1"] ==-1){
          tr[j,"A1"] = As[i]
        }else{
          tr[j,"A2"] = As[i]
        }
      }
    }
  }
  
  if(any(is.infinite(As))){return(list(Infinite=T))}
  
  
  for(i in 2:nrow(tr)){
    father=tr[i,2]
    j=which(tr[,1]==father)
    
    if(Rs[i]==-Inf){
      A=max(1e-300,(As[j]-As[i])/As[i])
      R=-log(exp((1/eta)*log(A))+1)
      R=max(R,Rmin)
      Rs[i]=R
    }
    int[i]=int[j]+Rs[i]
  }
  
  tr=cbind(tr,int)
  return(list(transform=tr,unsampled=unsampled,int=int,Rs=Rs,As=As,p_nUnsampled=p_nUnsampled,p_unsampled=p_unsampled,list_depth=list_depth,Infinite=F))
}

# change the number of unsampled nodes and the intervals for a given node set
change_int_eta=function(enhance,ind,tree,alpha,beta,eta,epsilon,lambdaN,aN,change=F,uns){
  N=tree$Nnode
  rm=c()
  rmR=c()
  unsampled=enhance$unsampled
  int=enhance$int
  p_nUnsampled=enhance$p_nUnsampled
  list_depth=enhance$list_depth
  continue=T
  if(is.null(ind)){
    rm=c()
  }else{
    n=extract.clade(tree,ind)$Nnode+1
    tr=enhance$transform
    node=which(tr[,1]==ind)
    if(unsampled[ind]>0){
      rm=(node-1):(node-unsampled[ind])
    }
    nodeR=node-length(rm)
    rmR=rm
    if(n==1){unsampled[ind]=-1
    }else if(ind==(N+2)){unsampled[ind]=-1}else{
      lambda=lambdaN[n]
      an=aN[n]
      p_nUnsampled=p_nUnsampled+(unsampled[ind]+1)*log(lambda+an)-(unsampled[ind])*log(lambda)
      u=runif(1)
      if(u<(0.4)){
        unsampled[ind]=unsampled[ind]-1-floor(rexp(1,0.2))
        ti=unsampled[ind]
      }else if(u>0.6){
        unsampled[ind]=unsampled[ind]+1+floor(rexp(1,0.2))
        ti=unsampled[ind]
      }else{
        if(unsampled[ind]>1){
          lrm=sample(min(unsampled[ind],maxChange),1)
          rmR=sample(rm,lrm)
          nodeR=node-lrm
          ti=lrm
        }else{
          ti=unsampled[ind]
        }
      }
      node=node-length(rm)
      if(unsampled[ind]<0){
        p_nUnsampled=-Inf
        continue=F
      }else{
        p_nUnsampled=p_nUnsampled-(unsampled[ind]+1)*log(lambda+an)+(unsampled[ind])*log(lambda)
        continue=T
      }
      
      if(is.na(unsampled[ind])) print(paste("a",alpha,"b",beta,"t",unsampled[ind],"lamb",lambda,"an",an,"sig",sigma,"n",n))
    }
  }
  if(! continue){
    enhance$p_nUnsampled=-Inf
    return(enhance)
  }
  
  tr=transform_eta(tree,unsampled,eta,list_depth = list_depth,change_depth = max(ind,0))
  p_unsampled=tr$p_unsampled
  list_depth=tr$list_depth
  tr=tr$M
  
  if(length(rm)>0){
    Rs=enhance$Rs[-rmR]
    depth=enhance$transform[-rm,"depth"]
    
  }else{
    Rs=enhance$Rs
    depth=enhance$transform[,"depth"]
    
  }
  if(! is.null(ind)){
    if(unsampled[ind]>-1){
      add=rep(-Inf,ti+1)
      if(nodeR<length(Rs)){
        Rs=c(Rs[1:(nodeR-1)],add,Rs[(nodeR+1):length(Rs)])
      }else{
        Rs=c(Rs[1:(nodeR-1)],add)}
      if(node<length(depth)){
        depth=c(depth[1:(node-1)],tr[node:(node+unsampled[ind]),"depth"],depth[(node+1):length(depth)])
      }else{
        depth=c(depth[1:(node-1)],tr[node:(node+unsampled[ind]),"depth"])
      }
    }
  }
  
  tr[,"depth"]=depth
  
  Rs[tr[,"father"]<(2*N+2)]=-Inf
  if(change){Rs[]=-1}
  Rs[1]=0
  int=c(0,rep(-Inf,nrow(tr)-1))
  As=rep(-1,nrow(tr))
  
  for(i in 1:N){
    do=(tr[,2]>=(2*N+2) & Rs==-Inf & tr[,4]==(i+1))
    samp=sample(10000,sum(do),replace=T)
    if(sum(do)>0){Rs[do]=uns[[i]][samp]}
  }
  
  nTr= nrow(tr)
  nodes=rep(T,nTr)
  while(sum(nodes)>0){
    i=which.min(depth)
    node=tr[i,"node"]
    if((tr[i,"A1"]==-1) | (tr[i,"A2"]==-1 & tr[i,"node"]<(2*N+2))){
      print("error in change int")
    }else{
      depth[i]=Inf
      if(tr[i,"A2"]>-1){
        As[i]=tr[i,"A1"]+tr[i,"A2"]
      }else{
        k=which(tr[,"father"]==node)
        As[i]=tr[i,"A1"]*(1+exp(eta*(log(1-exp(Rs[k]))-(Rs[k]))))
      }
      father=tr[i,"father"]
      nodes[i]=F
      if(father>-1) {
        j=which(tr[,"node"]==father)
        if(tr[j,"A1"] ==-1){
          tr[j,"A1"] = As[i]
        }else{
          tr[j,"A2"] = As[i]
        }
      }
    }
  }
  
  if(any(is.infinite(As))){return(list(Infinite=T))}
  
  
  for(i in 2:nrow(tr)){
    father=tr[i,2]
    j=which(tr[,1]==father)
    
    if(Rs[i]==-Inf){
      A=max(1e-300,(As[j]-As[i])/As[i])
      R=-log(exp((1/eta)*log(A))+1)
      R=max(R,Rmin)
      
      Rs[i]=R
    }
    int[i]=int[j]+Rs[i]
  }
  
  tr=cbind(tr,int)
  return(list(transform=tr,unsampled=unsampled,int=int,Rs=Rs,As=As,p_nUnsampled=p_nUnsampled,p_unsampled=p_unsampled,list_depth=list_depth,Infinite=F))
}

# and the inference function :

# run the inference of alpha and eta
# tree is the phylogeny with an additional field (tip.ab) with species abundances,
# epsilon is the minimum split size (as described i appendix 1),
# beta is the MLE of the beta statistic for the tree, 
# lambadN and aN are the vector given by the functions lambda_N and a_N, 
# niter is the number of iterations in the mcmc chain,  
# ini is a vector with the initial value of alpha and eta, V the initial scaling factors for the proposal functions,
# [ma,Ma] and [me,Me] are the intervals within which alpha and eta are allowed to be, 
# verbose is the number of iteration after which the state of the mcmc is printed (if silent == F)
# the proposal is scaled every Nadapt iterations (if at least NadaptMin iterations have been performed)
# chain is either NULL (and a new inference is being performed) or the output of mcmc_eta (and the chain is continued)
# An example of the use of this function is given in the end of this document
mcmc_eta=function(tree,epsilon,beta,ini=c(0,1),V=c(0.1,0.1),chain=NULL,
                  niter,verbose=10,silent=T,Nadapt=Inf,NadaptMin=10,
                  ma=-4,Ma=4,me=0.1,Me=10,proposal="bactrian",accOpt=0.3,lambdaN,aN){
  
  ND=rank(nodes.depths.ordonnes(tree))
  tree2=build.tree(ND)
  tree$edge.length=tree2$edge.length
  
  N=tree$Nnode
  sign=F
  
  proba.int=function(eta,enhanced){
    p=0
    if(any(enhanced$Rs)==-Inf){return(-Inf)}
    for(i in (N+2):(2*N+1)){
      ind=which(enhanced$transform[,"node"]==i)
      son=which(enhanced$transform[,"father"]==i)
      if(length(son)>0){
        R=enhanced$Rs[son[1]]
        Rm=log(1-exp(R))
        R=max(R,Rmin)
        Rm=max(Rm,Rmin)
        Ntip=enhanced$transform[ind,"ntip"]
        ntip=enhanced$transform[son[1],"ntip"]
        cond=max(1-exp(Ntip*R)-exp(Ntip*Rm),exp(Rmin))
        p=p+(ntip)*R+(Ntip-ntip)*Rm-log(cond)
      }
    }
    return(p)
  }
  
  proba.int2=function(eta,enhanced){
    p=0
    if(any(enhanced$Rs)==-Inf){return(-Inf)}
    for(i in (N+2):(2*N+1)){
      ind=which(enhanced$transform[,"node"]==i)
      son=which(enhanced$transform[,"father"]==i)
      if(length(son)>0){
        R=enhanced$Rs[son[1]]
        Rm=log(1-exp(R))
        R=max(R,Rmin)
        Rm=max(Rm,Rmin)
        Ntip=enhanced$transform[ind,"ntip"]
        ntip=enhanced$transform[son[1],"ntip"]
        p=p+(ntip+beta)*R+(Ntip-ntip+beta)*Rm
      }else{
        son=(tree$edge[tree$edge[,1]==i,2])
        A=max(tree$tip.ab[son[1]]/tree$tip.ab[son[2]],tree$tip.ab[son[2]]/tree$tip.ab[son[1]])
        R=1/(1+A^(1/eta))
        Rm=1-R
        R=max(log(R),Rmin)
        Rm=max(log(Rm),Rmin)
        p=p+(beta+1)*R+(1+beta)*Rm
      }
    }
    return(p)
  }
  
  posterior2=function(alpha,eta,enhanced){
    if(enhanced$Infinite | eta<0){return(-Inf)}
    Y=try(aux_lik(enhanced$transform,beta,alpha))
    if(inherits(Y,"try-error")){Y=-Inf}
    p=proba.int2(eta,enhanced)
    p2=enhanced$p_nUnsampled
    p3=enhanced$p_unsampled
    if(any(is.na(c(Y,p,p2,p3))) | any(is.nan(c(Y,p,p2,p3))) | any(is.infinite(c(Y,p,p2,p3)))){return(-Inf)}else{}
    return(Y+p-p3+p2)
  }
  
  posterior=function(alpha,eta,enhanced){
    if(enhanced$Infinite | eta<0){return(-Inf)}
    Y=try(aux_lik(enhanced$transform,beta,alpha))
    if(inherits(Y,"try-error")){Y=-Inf}
    p=proba.int(eta,enhanced)
    if(any(is.na(c(Y,p))) | any(is.nan(c(Y,p))) | any(is.infinite(c(Y,p)))){return(-Inf)}
    return(Y+p)
  }
  
  if (is.null(chain)){
    uns=lapply(1:(N),function(i){(-1*simulate.Yi(10000,epsilon,beta,i+1))})
    a=ini[1]
    e=ini[2]
    j=1
    post_actuel=-Inf
    int=enhance_eta(tree,a,beta,e,epsilon,lambdaN = lambdaN, aN=aN, uns=uns)
    post_actuel=posterior(a,e,int)
    while(post_actuel==-Inf){
      int=enhance_eta(tree,a,beta,e,epsilon,lambdaN = lambdaN, aN=aN, uns=uns)
      post_actuel=posterior(a,e,int)
    }
    Mcmc=c(a,e,sum(int$unsampled[(N+2):(2*N+1)]),post_actuel)
  }else{
    uns=chain$uns
    a=chain$a
    e=chain$e
    int=chain$int
    Mcmc=chain$mcmc
    j=nrow(chain$mcmc)
    V=chain$V
    post_actuel=chain$post_actuel
  }
  n=0
  for(i in 1:niter){
    if((i+j) %% Nadapt ==0 & (i+j)>NadaptMin){ 
      if(proposal=="uniform" | proposal=="bactrian"){
        acc1=sum(diff(Mcmc[-c(1:floor(nrow(Mcmc)-Nadapt)),1])!=0)/(Nadapt-1)
        if(acc1<0.001){
          V[1]= V[1]/100
        }else if(acc1>0.999){
          V[1]= V[1]*100
        }else{
          V[1] =  V[1] * (tan((pi/2)*acc1)/tan((pi/2)*accOpt))
        }
        acc2=sum(diff(Mcmc[-c(1:floor(nrow(Mcmc)-Nadapt)),2])!=0)/(Nadapt-1)
        if(acc2<0.001){
          V[2]= V[2]/100
        }else if(acc2>0.999){
          V[2]= V[2]*100
        }else{
          V[2] =  V[2] * (tan((pi/2)*acc2)/tan((pi/2)*accOpt))
        }
      }else{
        V[1]=max(1e-10,sqrt(var(Mcmc[-c(1:floor(nrow(Mcmc)/2)),1]))*2.36)
        V[2]=max(1e-10,sqrt(var(log(Mcmc[-c(1:floor(nrow(Mcmc)/2)),2])))*2.36)
      }
      if(any(is.na(V))) V=c(1,1)
    }
    
    #Random parameter drawing : 
    if(proposal=="uniform"){
      new_a = a + (runif(1,0,1)-0.5)*3.4641016*V[1]
      new_e =exp(log(e) + (runif(1,0,1)-0.5)*3.4641016*V[2])
      
    }else if (proposal=="bactrian"){
      new_a =a + rbactrian(1)*V[1]
      new_e =exp(log(e) + rbactrian(1)*V[2])
      
    }else{
      new_a=rnorm(1,a,V[1])
      new_e=rnorm(1,e,V[2])
    }
    
    new_int=int
    post_e=posterior2(a,e,new_int)
    if(post_e>-Inf){
      ks=(N+3):(2*N+1)
      for(k in ks){
        new_new_int=change_int_eta(new_int,k,tree,a,beta,e,epsilon,lambdaN = lambdaN, aN=aN, uns=uns)
        new_post=posterior2(a,e,new_new_int)
        
        if(new_post>=post_e){
          new_int=new_new_int
          post_e=new_post
        }else{
          u=runif(1,0,1)
          if (log(u)<(new_post-post_e)){
            new_int=new_new_int
            post_e=new_post
          }
        }}
      post_new_int=posterior(a,e,new_int)
      if(post_new_int>-Inf){
        int=new_int
        post_actuel=post_new_int
      }else{print("infinite")}
    }
    
    new_int=int
    if(new_a>Ma| new_a<ma){
      post=-Inf
    }else{
      post=posterior(new_a,e,int)
      
    }
    
    if(post>=post_actuel){
      a=new_a
      post_actuel=post
    }else{
      u=runif(1,0,1)
      if (log(u)<(post-post_actuel)){
        a=new_a
        post_actuel=post
      }
    }
    
    if(new_e>Me| new_e<me){
      post=-Inf
    }else{
      new_int=change_int_eta(int,NULL,tree,a,beta,new_e,epsilon,lambdaN = lambdaN, aN=aN, change=F, uns=uns)
      post=posterior(a,new_e,new_int)
    }
    if(post>=post_actuel){
      e=new_e
      int=new_int
      post_actuel=post
    }else{
      u=runif(1,0,1)
      if (log(u)<(post-post_actuel)){
        e=new_e
        int=new_int
        post_actuel=post
      }
    }
    
    Mcmc=rbind(Mcmc,c(a,e,sum(int$unsampled[(N+2):(2*N+1)]),post_actuel))
    
    n=n+1
    if(!silent){if(n==verbose){
      n=0
      print(paste(i,"alpha",Mcmc[i+j,1],"eta",Mcmc[i+j,2],"beta",beta,"post",Mcmc[i+j,4],"ntir",sum(int$unsampled[(N+2):(2*N+1)])))
    }}
  }
  colnames(Mcmc)=c("alpha","eta","nUnsampledled","log_post")
  return(list(mcmc=mcmc(Mcmc),post_actuel=post_actuel,a=a,int=int,tree=tree,V=V,e=e,uns=uns))
}


#### Example for the inference of alpha ####
run=F
if(run){
  seed=123
  ntip=30
  set.seed(seed)
  tree=simulate.tree(epsilon = 0.01,alpha = -1,beta = 0,N = ntip,equal.ab = T)
  beta=maxlik.betasplit(tree,up=10)$max_lik
  plot(tree)
  
  niter=1000
  lambdaN=lambda_N(epsilon = 0.01,beta = beta,N = ntip)
  aN=a_N(beta = beta,N = ntip)
  
  chain=mcmc_alpha(tree,epsilon=0.01,beta=beta,niter=niter,V = c(0.1),ini=c(0),
                   verbose = 100,silent = F,Nadapt = 200,NadaptMin = 200,
                   lambdaN=lambdaN,aN=aN)
  
  
  chain=mcmc_alpha(tree,epsilon=0.01,beta=beta,niter=20*niter,verbose = 500,silent = F,
                   chain = chain,Nadapt = 500,NadaptMin = 500,
                   lambdaN=lambdaN,aN=aN)
  
  thinned=mcmc(chain$mcmc[seq(5000,21000,10),])
  plot(thinned)
  da=density(thinned[,"alpha"])
  MPa=da$x[which.max(da$y)]
  print(MPa)
}

#### Example for the inference of alpha and eta ####

run=F
if(run){
  seed=123
  set.seed(seed)
  ntip=30
  tree=simulate.tree(epsilon = 0.001,alpha = -1,beta = 0,N = ntip,equal.ab = F,eta =1.5)
  beta=maxlik.betasplit(tree,up=10)$max_lik
  extinctions = rank(tree$tip.ab)
  tree$tip.label = rep(".", length(tree$tip.label))
  plot.phylo(tree, show.node.label=TRUE, cex=order(extinctions, seq(1,(tree$Nnode+1)))/((tree$Nnode+1)/6), adj=0.1)  
  
  niter=1000
  lambdaN=lambda_N(epsilon = 0.001,beta = beta,N = ntip)
  aN=a_N(beta = beta,N = ntip)
  
  chain=mcmc_eta(tree,epsilon=0.001,beta=beta,V = c(0.1,0.1),niter=niter,ini=c(0,1),
                 verbose = 5,silent = F,Nadapt = 200,NadaptMin = 200,
                 lambdaN=lambdaN,aN=aN)
  #The initialisation of the mcmc is quiet long because we begin by drawing many unsampled intervals. When this is done it
  #should get quicker. 
  
  chain=mcmc_eta(tree,epsilon=0.001,beta=beta,niter=20*niter,verbose = 20,silent = F,
                 chain = chain,Nadapt = 500,NadaptMin = 500,
                 lambdaN=lambdaN,aN=aN)
  
  thinned=mcmc(chain$mcmc[seq(5000,21000,10),])
  plot(thinned)
  da=density(thinned[,"alpha"])
  MPa=da$x[which.max(da$y)]
  de=density(log(thinned[,"eta"]))
  MPe=exp(de$x[which.max(de$y)])
  print(MPa)
  print(MPe)
}

chain=mcmc_eta(tree,epsilon=0.001,beta=beta,niter=20*niter,verbose = 20,silent = F,
               chain = chain,Nadapt = 500,NadaptMin = 500,
               lambdaN=lambdaN,aN=aN)
