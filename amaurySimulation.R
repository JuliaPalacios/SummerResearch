# This code runs the model presented in "Ranked tree shapes, 
# non-random extinctions and the loss of phylogenetic diversity",
# by Odile Maliet, Fanny Gascuel and Amaury Lambert (submitted
# to Systematic Biology). Its allows to simulate random ranked tree shapes
# endowed with random numbers summing at one at the tips (interpreted as
# relative species abundances), as a function of parameters beta, alpha and eta
# (as detailed in the article and Supplementary Material).
# The code uses the algorithm presented in the Appendix 1, available in the online
# Supplementary Material.

######functions############
library(cubature)
library(ape)
library(pbapply)
library(picante)
library(apTreeshape)

# Simulates the parameter lambda.epsilon of the online appendix
# epsilon and beta are the parameters of the online appendix
# n is the number of marks in the interval
lambda.epsilon=function(epsilon,beta,n,upper=10000){
  f=function(x){
    sapply(x,function(e){exp(-(beta+n+1)*e)*(1-exp(-e))^beta})
  }
  # integrate(f,epsilon,upper)
  int=2*adaptIntegrate(f, lowerLimit = epsilon, upperLimit = upper)$integral
  while(upper>0.1 & int==0){
    upper=upper/1.5
    int=2*adaptIntegrate(f, lowerLimit = epsilon, upperLimit = upper)$integral
  }
  return(int)
}

lambda_N=function(epsilon,beta,N,upper=10000){
  c(0,sapply(2:N, function(x){lambda.epsilon(epsilon,beta,x,upper)}))
}

# Simulates the random variables Y_i (see online appendix)
# N is the number of simulated Y_i
# epsilon and beta are the parameters of the online appendix
# n is the number of marks in the interval
simulate.Yi=function(N,epsilon,beta,n){
  res=rep(0,N)
  if(N<1) return(c())
  
  if(beta<=0){
    g=function(x){((1-exp(-epsilon))^beta)*exp(-(beta+n+1)*epsilon)}  
    
    f=function(x){
      sapply(x,function(e){(1-exp(-e))^beta})
    }
    
    fun=function(i){
      continue=T
      while(continue){
        e=rexp(1,(beta+n+1))+epsilon
        u=runif(1,0,1)
        if(u<(f(e)/g(e))){
          return(e)
          continue=F
        }
      }
    }
    
    res=sapply(1:N,fun)
    
  }else{
    
    simulate.g=function(N,beta,n){
      res=rep(0,N)
      A=(beta+n+1)
      M=((beta/(A+beta))^beta)*(1+beta/A)^(-A)
      for(i in 1:N){
        u=runif(1,0,1)
        if(u<(-log(M)/(1-log(M)))){
          res[i]=runif(1,0,-log(M)/A)
        }else{
          res[i]=-log(M)/A+rexp(1,A)
        }
      }
      return(res)
    }
    
    f=function(x){
      A=beta+n+1
      M=((beta/(A+beta))^beta)*(1+beta/A)^(-A)
      sapply(x,function(e){exp(-A*e)*(1-exp(-e))^beta})
    }
    
    g=function(x){
      A=beta+n+1
      M=((beta/(A+beta))^beta)*(1+beta/A)^(-A)
      sapply(x,function(y){
        min(exp(-A*y),M)
      })
    }
    
    fun=function(i){
      continue=T
      while(continue){
        e=simulate.g(1,beta,n)
        if(e>epsilon){
          u=runif(1,0,1)
          if(u<(f(e)/g(e))){
            return(e)
            continue=F
          }}
      }
    }
    res=sapply(1:N,fun)
  }
  
  return(res)
}

# Simulates the time before the splitting of the marks into two different fragments,
# and the size of the fragment at this time.
# epsilon, alpha, beta and eta are parameters introduced in the article and the online Appendix 1
# x is the initial size of the interval, n is the initial number of marks in the interval
# ab control whether relative abundances are computed or not (if one wants to use a simple field-of-bullets model of extinctions)
# x.ab is the initial abundance
simulate.Tau.X=function(epsilon,x,alpha,beta,n,ab=F,eta=1,x.ab=1,lambda=NULL){
  if(is.null(lambda)){
    lambda=lambda.epsilon(epsilon,beta,n)
  }else{
    lambda=lambda[n]
  }
  an=sum(sapply(1:(n-1),function(i){beta(beta+1+i,beta+1+n-i)/((n+1)*beta(1+i,1+n-i))}))
  sigma=rexp(1,an)
  N=rpois(1,lambda*sigma)
  Y=c(simulate.Yi(N,epsilon,beta,n))
  Z=cumsum(c(0,Y))
  si=diff(c(0,sort(runif(N,0,sigma)),sigma))
  Tau=x^(-alpha)*(sum(exp(alpha*Z)*si))
  X=x*exp(-Z[length(Z)])
  if (ab){
    if(N==0){
      X.ab=x.ab
      return(c(Tau,X,X.ab))
    }
    
    if(eta==0){
      X.ab=x.ab/(2^N)
    }else{
      Y2=sum(log(exp(-Y*eta)+(1-exp(-Y))^eta))
      X.ab=x.ab*exp(-Z[length(Z)]*eta)/exp(Y2)
    }
    return(c(Tau,X,X.ab))
  }
  return(c(Tau,X))
}

# Simulates the number of marks in the interval and the size of the interval after splitting
# beta is a parameter of the model
# n is the initial number of marks in the interval
simulate.R.K=function(beta,n){
  prob=sapply(1:(n-1),function(i){beta(beta+1+i,beta+1+n-i)/beta(1+i,1+n-i)})
  k=sample(1:(n-1),size=1,prob=prob)
  r=rbeta(1,beta+1+k,beta+1+n-k)
  return(c(r,k))
}

# Simulates the size of the interval after splitting, knowing the number of marks in this interval
# beta is a parameter of the model
# n is the initial number of marks in the interval
# K is the number of marks in the left subinterval after splitting
simulate.R=function(beta,n,K){
  r=rbeta(1,beta+1+K,beta+1+n-K)
  return(c(r,K))
}

# Binds two trees a1 and a2, each of them being respectively at a distance root1 and root2 from the root
# ab controls whether relative abundances are recorded in the tree
bind.trees=function(a1,a2,root1,root2,ab=F){
  if(a1$Nnode==0){return(a2)
  }else if (a2$Nnode==0){return(a1)
  }else{
    a1$root.edge=root1
    a2$root.edge=root2
    a=a1+a2
    if(ab){a$tip.ab=c(a1$tip.ab,a2$tip.ab)}
    return(a)
  }
}

# Simulates the ranked topology 
# epsilon, alpha, beta and eta are parameters of the model, as introduced in the article and online Appendix 1
# N is the number of tips of the tree
# if equal.ab is set to TRUE, all the tips have equal abundances. If not, the model simulates relative abundances 
simulate.tree=function(epsilon,alpha,beta,N,equal.ab=T,eta=1,lambda=NULL){
  aux=function(x,n,x.ab=1){
    if(n>1){
      v1=simulate.R.K(beta,n)
      n1=v1[2]
      n2=n-n1
      x1=v1[1]*x
      x2=x-x1
      if(n1>1){
        if(equal.ab){v2.1=simulate.Tau.X(epsilon,x1,alpha,beta,n1,lambda=lambda)
        }else{v2.1=simulate.Tau.X(epsilon,x1,alpha,beta,n1,ab=T,eta = eta,x.ab=(x1^eta)/(x1^eta+x2^eta),lambda=lambda)}
        tree1=aux(v2.1[2],n1)
        if(!equal.ab){
          tree1$tip.ab=v2.1[3]*tree1$tip.ab
        }
        root1=v2.1[1]
      }else{
        tree1=list(edge=matrix(c(2,1),1,2),edge.length=0,Nnode=1,tip.label=NA,tip.ab=c(1))
        if(!equal.ab){
          tree1$tip.ab=((x1^eta)/(x1^eta+x2^eta))
        }
        class(tree1)="phylo"
        
        root1=1
      }
      if(n2>1){
        if(equal.ab){v2.2=simulate.Tau.X(epsilon,x2,alpha,beta,n2,lambda=lambda)
        }else{v2.2=simulate.Tau.X(epsilon,x2,alpha,beta,n2,ab=T,eta = eta,x.ab=(x2^eta)/(x1^eta+x2^eta),lambda=lambda)}
        tree2=aux(v2.2[2],n2)
        if(!equal.ab){
          tree2$tip.ab=v2.2[3]*tree2$tip.ab
        }
        root2=v2.2[1]
      }else{
        tree2=list(edge=matrix(c(2,1),1,2),edge.length=0,Nnode=1,tip.label=NA,tip.ab=c(1))
        if(!equal.ab){
          tree2$tip.ab=((x2^eta)/(x1^eta+x2^eta))
        }
        class(tree2)="phylo"
        
        root2=1
      }
      a=bind.trees(tree1,tree2,root1,root2,ab=T)
      A=collapse.singles(a)
      A$tip.ab=a$tip.ab
      return(A)
    }else{
      tree=list(edge=matrix(c(2,1),1,2),edge.length=1,Nnode=1,tip.label=NA,tip.ab=1)
      class(tree)="phylo"
      return(tree)
    }
  }
  
  tree=aux(1,N)
  order=c(rep(N,N),rank(node.depth.edgelength(tree)[(N+1):(2*N-1)]))
  for(i in 1:(2*N-2)){
    tree$edge.length[i]=order[tree$edge[i,2]]-order[tree$edge[i,1]]
  }
  return(tree)
}

myFavoriteTree = simulate.tree(2,100,20,10,equal.ab=FALSE,eta=1,lambda=NULL)
ggtree(myFavoriteTree)
# Builds the tree using node depths as in the birth-??death model and coalescent point process;
# "The shape and probability of reconstructed phylogenies",
# by Amaury Lambert and Tanja Stadler (Theoretical Population Biology, 2013)
# H is a vector of node depths (the order matters)
# tip.lab is a vector of tip names
build.tree=function(H,tip.lab=rep(NA,length(H)+1)){
  if (length(H)==0) {
    tree=list(edge=matrix(c(2,1),1,2),edge.length=0,Nnode=1,tip.label=tip.lab)
    class(tree)="phylo"
    return(tree)
  }else{
    H.crown=max(H)
    split=which.max(H)
    t1=tip.lab[1:split]
    t2=tip.lab[(split+1):length(tip.lab)]
    if (split>1) {H1=H[1:(split-1)]
    }else{H1=c()}
    if (split<length(H)){H2=H[(split+1):length(H)]
    }else {H2=c()}
    
    a1=build.tree(H1,t1)
    a2=build.tree(H2,t2)
    root1=H.crown-max(c(H1,0))
    root2=H.crown-max(c(H2,0))
    return(collapse.singles(bind.trees(a1,a2,root1,root2)))
  }
}



# Gives the node depths of the tree, in the order of node labels 
depth=function(tree){
  d=node.depth.edgelength(tree)
  d[1]-d[(tree$Nnode+2):(2*tree$Nnode+1)]
}

# Gives the node depths of the tree, in the order needed to reconstruct the tree using build.tree
nodes.depths.ordonnes=function(tree){
  aux=function(phylo,noeud=phylo$Nnode+2){
    if(phylo$Nnode<=1){return(max(phylo$edge.length))
    }else{
      fils=phylo$edge[phylo$edge[,1]==noeud,][,2]
      if (fils[1]<(phylo$Nnode+2)){
        ad=extract.clade(phylo,fils[2])
        return(cbind(max(phylo$edge.length[phylo$edge[,1]==noeud]),aux(ad)))
      }else if (fils[2]<(phylo$Nnode+2)){
        ag=extract.clade(phylo,fils[1])
        return(cbind(aux(ag),max(phylo$edge.length[phylo$edge[,1]==noeud])))
      }else{
        ag=extract.clade(phylo,fils[1])
        ad=extract.clade(phylo,fils[2])
        hg=aux(ag)
        hd=aux(ad)
        h=max(hg)+phylo$edge.length[phylo$edge[,2]==fils[1]]
        return(cbind(hg,h,hd))
      }
    }
  }
  H=aux(tree)
  return(H)
}

# Simulates the ranked topology with node depths simulated using the Kingman's coalescent model
# epsilon, alpha, beta and eta are parameters of the model, as introduced in the article and online Appendix 1
# N is the number of tips of the tree
# if equal.ab is set to TRUE, all the tips have equal relative abundances. If not, relative abundances are simulated with the model
simulate.kingman=function(epsilon,alpha,beta,N,equal.ab=T,eta=1,lambda=NULL){
  tree=simulate.tree(epsilon,alpha,beta,N,equal.ab,eta,lambda=lambda)
  ab=tree$tip.ab
  order=rank(nodes.depths.ordonnes(tree))
  depths=cumsum(rexp(N-1,sapply((N:2),function(i){i*(i-1)/2})))
  tree=build.tree(depths[order],1:N)
  tree$tip.ab=ab
  return(tree)
}

# Orders the tips which will sequentially get extinct/be sampled, based on relative species abundances,
# as explained in the main text of the article.
# First species in output = first to go extinct
# N is the number of tips
# tree is the phylogeny
# If equal.ab is set to TRUE, all species have the same probability to go extinct first;
# if not, the order of species extinctions is determined deterministically in the increasing order of relative abundances at the tips.
get.extinction.list=function(N,tree,equal.ab=T){
  
  tip.ab=tree$tip.ab
  names(tip.ab) = tree$tip.label
  tips.order = NULL
  if (equal.ab == T) {
    tips.order = sample(names(tip.ab), N)
    
  } else {  
    ab_null = c(which(is.na(tip.ab)), which(tip.ab==0))
    if (length(ab_null)>0) {
      ab_null_order = sample(ab_null, length(ab_null))
      tip_not_null = tip.ab[-ab_null]
      tips.order = c(names(ab_null_order), names(tip_not_null)[order(tip_not_null,names(tip_not_null))]) 
    } else {
      tips.order = names(tip.ab)[order(tip.ab,names(tip.ab))] 
    }
  } 
  return(tips.order)
}

# Simulates node depths in the birth-death model, conditionned by the number
# of tips and the age of the root (using the expression in "The shape and 
# probability of reconstructed phylogenies", by Amaury Lambert and Tanja Stadler 
# (Theoretical Population Biology, 2013))
# N is the number of tips, b the birth rate, d the death rate, tmax the maximum age of the root
yule.lengths=function(N,b,d,tmax){
  r=b-d
  
  f=function(u){
    if (r==0){
      u/(b*(1-u))
    }else{
      log((r*u/(b*(1-u)))+1)/r
    }
  }
  
  X=c(tmax)
  n=1
  while (n<N){
    t=f(runif(1,0,1))
    if (t<tmax){X=c(X,t);n=n+1}
  }
  
  return(X[2:N])
}

# Simulates a phylogeny in a birth-death model, conditionned by the number
# of tips and the maximum age of the root 
# N is the number of tips, b the birth rate, d the death rate, tmax the maximum
# age of the root
yule=function(N,b,d,tmax){
  build.tree(yule.lengths(N,b,d,tmax))
}

myTree = build.tree(yule.lengths(5,2,2,3))
ggtree(myTree)

# Simulates the ranked topology, with node depths as in the birth-death model
# epsilon, alpha, beta and eta are parameters of the model, as introduced in the article and online Appendix 1
# N is the number of tips of the tree
# b the birth rate, d the death rate, tmax the age of the root
# if equal.ab is set to TRUE, all the tips have equal relative abundances. If not, relative abundances are simulated with the model
simulate.yule=function(epsilon,alpha,beta,N,b,d,tmax=Inf,equal.ab=T,eta=1,lambda=NULL){
  tree=simulate.tree(epsilon,alpha,beta,N,equal.ab,eta,lambda=lambda)
  ab=tree$tip.ab
  order=rank(nodes.depths.ordonnes(tree))
  depths=sort(yule.lengths(N,b,d,tmax))
  tree=build.tree(depths[order],1:N)
  tree$tip.ab=ab
  return(tree)
}

# Computes the proportion of conserved phylogenetic diversity as a function of the proportion
# of conserved species, in trees simulated by the model
# epsilon, alpha, beta and eta are parameters of the model, as introduced in the article and online Appendix 1
# ntree is the number of trees simulated
# N is the number of tips of the trees
# sample.frac is the fractions of tips for which we want to compute the conserved PD 
# If lengths="yule", node depths correspond to a birth-death process (with birth rate b, death rate d, root age tmax),
# if lengths="kingman", node depths correspond to the Kingman's coalescent model 
# If equal.ab is set to TRUE, all species have the same probability to go extinct first, 
# if not, the order of species extinctions is determined deterministically in the order of increasing relative abundances.
get.PD.sample=function(epsilon,beta,alpha,N,sampl.frac,ntree,equal.ab,eta,lengths="yule",b=1,d=0){
  lambda=lambda_N(epsilon,beta,N,upper=10000)
  PD=pbsapply(1:ntree,function(i){
    if(lengths=="yule"){
      tree=simulate.yule(epsilon,alpha,beta,N,b = b,d=d,equal.ab = equal.ab,eta = eta,lambda=lambda)
    }else if(lengths=="kingman"){
      tree=simulate.kingman(epsilon,alpha,beta,N,equal.ab = equal.ab,eta = eta,lambda=lambda)
    }
    M=max(depth(tree))
    tip=get.extinction.list(N,tree,equal.ab)
    PD.ech=sapply(1:length(sampl.frac),function(j){
      if (round(sampl.frac[j]*N)==N){sum(depth(tree))+M
      } else if (round(sampl.frac[j]*N) == 0) {
        0
      } else if (round(sampl.frac[j]*N) == 1) {
        M
      } else {
        nd=depth(drop.tip(tree,tip[1:(N-round(sampl.frac[j]*N))]))
        sum(nd)+M}})
    return(PD.ech/(sum(depth(tree))+M))
  })
  rownames(PD) = sampl.frac
  colnames(PD) = seq(1,ntree)
  return(PD)
}

# Computes the maximum likelihood estimate of the parameter beta in trees simulated by the model, 
# as a function of the proportion of conserved species 
# epsilon, alpha, beta and eta are parameters of the model, as introduced in the article and online Appendix 1
# ntree is the number of trees simulated
# N is the number of tips of the trees
# sample.frac is the fractions of tips for which we want to compute the conserved PD 
# If lengths="yule", node depths correspond to a birth-death process (with birth rate b, death rate d, root age tmax),
# if lengths="kingman", node depths correspond to the Kingman's coalescent model 
# If equal.ab is set to TRUE, all species have the same probability to go extinct first, 
# if not, the order of species extinctions is determined deterministically in the order of increasing relative abundances.
get.tree.beta=function(epsilon,beta,alpha,N,sampl.frac,ntree,equal.ab,eta,lengths="yule",b=1,d=0){
  beta=pbsapply(1:ntree,function(i){
    if(lengths=="yule"){
      tree=simulate.yule(epsilon, alpha, beta, N, b=b, d=d, equal.ab=equal.ab, eta=eta)
    }else if(lengths=="kingman"){
      tree=simulate.kingman(epsilon, alpha, beta, N, equal.ab=equal.ab, eta=eta)
    }
    tip=get.extinction.list(N,tree,equal.ab)
    return(sapply(1:length(sampl.frac),function(j){
      if (round(sampl.frac[j]*N)==N){maxlik.betasplit(as.treeshape(tree), up=10, remove.outgroup=F, confidence.interval="none")$max_lik
      } else if (round(sampl.frac[j]*N)==0 || round(sampl.frac[j]*N)==1) {
        NA
      } else {
        maxlik.betasplit(as.treeshape(drop.tip(tree,tip[1:(N-round(sampl.frac[j]*N))])), up=10, remove.outgroup=F, confidence.interval="none")$max_lik       
      }
    }))
  })
  rownames(beta) = sampl.frac
  colnames(beta) = seq(1,ntree)
  return(beta)
}


##### Example ####
run=F

if(run){
  # Exemple of results on the loss of phylogenetic diversity (1 - proportion of conserved PD)
  # as a function of the proportion of extinction species (1- proportion of conserved species), in trees simulated by the model
  sampl.frac=seq(0,1,by=0.01) # fractions of tips for which we want to compute the conserved PD 
  beta = 0   # tree balance (ranked tree shape parameter)
  alpha = 0  # correlation clade size-age (ranked tree shape parameter)
  eta = 1    # correlation clade size-frequency (distribution of species extinction risks)
  N = 100    # number of tips of the trees
  ntree = 100 # number of simulation replicates
  # we consider here node depths as in the birth-death process, with b=1 and d=0 (pure-birth)
  # epsilon is set to 0.01 ; it must be small enough, but smaller epsilon does not affect the results
  
  PD = get.PD.sample(epsilon=0.01,beta=beta,alpha=alpha,N=N,sampl.frac=sampl.frac,ntree=ntree,equal.ab=F,eta=eta,lengths="yule",b=1,d=0)
  probs = c(0.025, 0.975, 0.5)  
  PD_stats = as.vector(t(sapply(1:nrow(PD), function(i){quantile(PD[i,], probs)}))) 
  par(mgp=c(2.2, 0.8, 0))
  par(mar=c(4, 4, 1, 1))
  plot(1, type="n", xlab="Fraction of extinct species, p", ylab="PD loss", ylim=c(0,1), xlim=c(0,1))
  points(c(0,1),c(0,1),t="l",col=grey(0),lty=3)
  # plot 95% confidence intervals (grey area)
  polygon(c(1-sampl.frac, rev(1-sampl.frac)), c(1-PD_stats[(1:length(sampl.frac))], rev(1-PD_stats[((length(sampl.frac)+1):(2*length(sampl.frac)))])), border=NA, col=grey(0.7))
  # plot median value (black line) 
  points(1-sampl.frac, 1-PD_stats[((2*length(sampl.frac)+1):(3*length(sampl.frac)))],t="l")
  
}