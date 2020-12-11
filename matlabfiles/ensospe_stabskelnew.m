function [wr,wi,vects]=ensospe_stabskelnew(kk,ekr,QQ,GG,HH,soref)
%
% Skeleton model (deterministic or stochastic)
% by Sulian Thual
% 
% computes stability on K,Rm,Zm,Am (defined on Hermite functions)
% Adim is the one from MS2011

% Input:
% - kk=desired wavenumber, kk=j*(2pi/L) with j integer (adim)
% - parameters...
% - nyk=1 here, numbers of Hermites (for Am,Zm), from m=0
% - nrm=1  here, number of Rossbys, from m=1
% - Here it is assumed so(x,y)=so(y)  i.e. homogeneous zonally
% - 
% Output:
% [wr,wi,Vects]:
% - eigenvalues(nmods)
% - vects(nxs,nmods)
%
% Note: assumes identical dissipation ekr on all variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RCE state
MMaeq=soref/HH/0.7511;

%
% Compute Matrix for : dtX + AA X + CC dxX=0
% X=(K, Rm, Am, Zm ) with all hermite-strips concatenation   
% Rm starts from m=1
% Am, Zm starts from m=0
nxs=4; % size of X
ii=complex(0,1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute CC:
CC=zeros(nxs,nxs);
% Compute AA:
AA=zeros(nxs,nxs);
%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Kelvin: (dt+e)K+dxK=-Ha/2
AA(1,1)=ekr;%K:K
CC(1,1)=1;% K
AA(1,4)=HH/2;%K:A
%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Rossby: (dt+e)R-dxR/3=-Ha/3
AA(2,2)=ekr;%R:R
CC(2,2)=-1/3;% R
AA(2,4)=HH/3;%R:A
%
%%%%%%%%%%%%%%%%%%%%%%%%%
% Z: dtZ -(QQ-1)*HH*A=0
AA(3,4)=(1-QQ)*HH;%Zm:Am (Ok)
AA(3,3)=ekr;

%%%%%%%%%%%%%%%%%%%%%%%%%
% A: dtA=GG*MMA*Q
% Q=Z-QQ*theta=Z+QQ*K+QQ*R
% dtA=GG*MMA*Z+GG*MMA
AA(4,1)=-GG*MMaeq*QQ;
AA(4,2)=-GG*MMaeq*QQ;
AA(4,3)=-GG*MMaeq;
AA(4,4)=ekr;

%
%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%
% Compute stability
% System: dtX + AA X + CC dxX=0
% Solution : X=Xo*exp(i(kx-wt))
% -iw*X + AA*X + ik*CC*X=0
% dtX = (-AA -i*k*CC)X=MX, eigenvalues l(k)
% w(k)=i*l(k)
MM=AA-ii*kk*CC;
[V,d]=eig(MM);
term=size(d);
nmod=term(1); % can have less modes than nxs
vects=zeros(nxs,nmod);
eigenvals=zeros(nmod,1);
for ixs=1:nmod; 
eigenvals(ixs,1)=d(ixs,ixs); 
end
vects=V;
wr=real(ii*eigenvals);% w(k)=i*l(k)
wi=imag(ii*eigenvals);% w(k)=i*l(k)

% ORDER BY DECREASING SPEED:
%1=Kelvin dry (50 ms-1)
%2=MJO mode (5 ms-1)
%3=Rossby moist mode (-5 ms-1)
%4=Rossby dry mode (-15 ms-1)
if 1==1
cr=wr/kk;
[term,index]=sort(cr,'descend');% sort by decreasing speed
wr=wr(index);
wi=wi(index);
vects(:,:)=vects(:,index);
end

% Notes to recall
%wi=imag(ii*eigenvals);
%cr=wr/kk;
% Solution : X=Xo*exp(i(kx-wt))
%cos(2pi*t/T) has period T
% frequency=2pi/T=wr 
%Cycle per day=1/T=wr/2pi
% wavenumber=2pi/L=k=j*(2pi/40000km)
%wr=wr/ta*oneday;
%wi=wi/ta*oneday;

% Parameters (adim may be wrong !)
%onehour=3600.;% one hour in seconds
%oneday=onehour*24.;% one day in seconds
%onekm=1000.;% one km in meters
%xa=1500.*onekm; %x adim (meters), cf BielloMajda2004 
%dx=625.*onekm/xa;% dx grid step (adim)
%ta=8*onehour; %t adim (seconds)














