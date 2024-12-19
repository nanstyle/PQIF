install.packages("fda")
library("MASS");library("fda")
f=function(t){
  t=as.matrix(t)
  f=t*(as.numeric(t>0))
  return(f)
}
bp=function(t,kn,k){
  basis=bsplineS(t,kn,k,0)
  return(basis)
}
dbp=function(t,kn,k){
  basis=bsplineS(t,kn,k,1)
  return(basis)
}
plam=function(lam,theta){
  ID1=as.numeric(theta<=lam)
  ID2=as.numeric(theta>lam)
  plam=lam*(ID1+(f(3.7*lam-theta)/(2.7*lam))*ID2)
  return(plam)
}
pqif=function(y,x,z,N,T,l,Nknot,maxit,lam,derta0){
 M=N*T
 RR=sqrt(N)
 w=matrix(0,nrow=N,ncol=N)
  for(i1 in 1:(N-1)){
    if(i1%%RR!=0){
      w[i1,(i1+1)]=1;w[(i1+1),i1]=1}
  }
  for(i2 in 1:(N-RR)){
    w[i2,(i2+RR)]=1;w[(i2+RR),i2]=1
  }
  for(i3 in 1:N) {
    if(sum(w[i3, ])>0){
      w[i3, ]=w[i3, ]/sum(w[i3, ])}
  }
 W=kronecker(w,diag(1,T))
 IT<-t(rep(1,T));ET <- diag(1,T)-(t(IT)%*%IT)/T
 PNT <- kronecker(diag(1,N),ET)
 ywan <- PNT%*%y;xwan <- PNT%*%x
 nz=ncol(z)
 theta=matrix(rep(1,nz),nz)/sqrt(nz)
 tn=z %*% theta
 kn=c(min(tn),quantile(tn,(1:Nknot)/(1+Nknot)),max(tn))
 b=bp(tn,kn,l);db=dbp(tn,kn,l)
 bwan <- PNT%*%b;dbwan <- PNT%*%db
 derta=derta0
 Y=ywan-derta*W%*%ywan
 wxwan=W%*%xwan
 X=xwan-derta*wxwan
 beta=ginv(t(X)%*%X)%*%t(X)%*%Y
 wbwan=W%*%bwan
 B=bwan-derta*wbwan
 gamma=ginv(t(B)%*%B)%*%t(B)%*%(Y-X%*%beta)
 wdbwan=W%*%dbwan
 para=rbind(derta,beta,theta,gamma);ng=nrow(gamma);np=nrow(para);nx=ncol(X)
 I1=diag(1,T,T)
 I2=diag(0,T,T)
  for(k in 1:T-1){
    I2[k,k+1]=1
    I2[k+1,k]=1}
 I3=diag(0,T,T)
   I3[1,1]=1
   I3[T,T]=1
 L1=matrix(0,T-1,T)
  for(t in 1:(T-1)){
    L1[t,t]=-1;L1[t,t+1]=1}
 L2=ginv(L1)
 for(iter in 1:maxit){
  mu=X%*%beta+B%*%gamma
  s=0
  for(m in 1:M){
   s=s+(Y[m]-mu[m])^2/(M-np)}
  cn=0;gn=0;gndot=0
  for(i in 1:N){
    yi=Y[((i-1)*T+1):(i*T),]
    zi=z[((i-1)*T+1):(i*T),]
    bi=bwan[((i-1)*T+1):(i*T),]
    xi=xwan[((i-1)*T+1):(i*T),]
    tni=tn[((i-1)*T+1):(i*T),]
    mui=mu[((i-1)*T+1):(i*T),]
    wxi=wxwan[((i-1)*T+1):(i*T),]
    wbi=wbwan[((i-1)*T+1):(i*T),]
    dbi=dbwan[((i-1)*T+1):(i*T),]
    wdbi=wdbwan[((i-1)*T+1):(i*T),]
    muidot=cbind(-wxi%*%beta-wbi%*%gamma,xi-derta*wxi,
                 diag(as.vector(dbi%*%gamma))%*%zi-derta*diag(as.vector(wdbi%*%gamma))%*%zi,
                 bi-derta*wbi)
    ai=s*diag(1,T,T)
    gi1=t(L1%*%muidot)%*%t(L2)%*%sqrt(ginv(ai))%*%I1%*%sqrt(ginv(ai))%*%L2%*%L1%*%(yi-mui)
    gi2=t(L1%*%muidot)%*%t(L2)%*%sqrt(ginv(ai))%*%I2%*%sqrt(ginv(ai))%*%L2%*%L1%*%(yi-mui)
    gi3=t(L1%*%muidot)%*%t(L2)%*%sqrt(ginv(ai))%*%I3%*%sqrt(ginv(ai))%*%L2%*%L1%*%(yi-mui)
    gi=rbind(gi1,gi2,gi3)
    gidot1=-t(L1%*%(muidot))%*%t(L2)%*%sqrt(ginv(ai))%*%I1%*%sqrt(ginv(ai))%*%L2%*%
           (L1%*%muidot)
    gidot2=-t(L1%*%(muidot))%*%t(L2)%*%sqrt(ginv(ai))%*%I2%*%sqrt(ginv(ai))%*%L2%*%
           (L1%*%muidot)
    gidot3=-t(L1%*%(muidot))%*%t(L2)%*%sqrt(ginv(ai))%*%I3%*%sqrt(ginv(ai))%*%L2%*%
           (L1%*%muidot)
    gidot=rbind(gidot1,gidot2,gidot3)
    gn=gn+gi/N
    cn=cn+gi%*%t(gi)/N
    gndot=gndot+gidot/N
 }
   H=c(rep(0,np))
   for(i7 in 1:np){
    H[i7]=plam(lam,abs(para[i7]))/(abs(para[i7]))
   }
 PH=diag(H)
 Qndot=t(gndot)%*%ginv(cn)%*%gn+PH%*%para
 Qn2dot=(t(gndot)%*%ginv(cn)%*%gndot)+PH
 para=para-ginv(Qn2dot)%*%Qndot
 derta=para[1]
 maxvalue=1;minvalue=-1;
 derta=min(maxvalue, max(minvalue,derta))
 beta=matrix(para[2:(nx+1)],nx)
 theta1=para[(nx+2):(1+nx+nz)]
 theta=matrix(theta1/sqrt(sum(theta1^2)),nz)
 gamma=matrix(para[(nx+nz+2):np],np-nx-nz-1)
 if(sum(abs(ginv(Qn2dot)%*%Qndot)) <= 1e-05)
 break
 tn=z %*% theta
 kn=c(min(tn),quantile(tn,(1:Nknot)/(1+Nknot)),max(tn))
 b=bp(tn,kn,l)
 db=dbp(tn,kn,l)
 bwan <- PNT%*%b
 dbwan <- PNT%*%db
 wbwan=W%*%bwan
 wdbwan=W%*%dbwan
 B=bwan-derta*wbwan
 X=xwan-derta*wxwan
 Y=ywan-derta*W%*%ywan
 Qn=t(gn)%*%ginv(cn)%*%gn
 df=sum(diag(ginv(Qn2dot)%*%(Qn2dot-PH)))
}
return(rbind(derta,beta,theta,gamma,lam))
}