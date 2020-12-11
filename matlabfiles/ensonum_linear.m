% by Sulian Thual

% compute linear matrix of the ENSO-MJO model

% Notes:
% -This is only used for initial conditions at the beginning of a simulation
% -The method of lines is used here in both the ocean, intraseasonal atmosphere and interannual atmosphere
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model Matrix 
% Compute Matrix for : dtX  = MM X 
% X={Ka,Ra,A,Z,Ko,Ro,T} 
% in that order (projected on the ocean/atm parabolic cylinder functions).
% we use variable change:
%z=q+QQ*o=q-QQ*(Ka+Ra)
%dtz=(1-Q)*(-Ha)
%dta=(GG*Ma)*q=(GG*Ma)*(z+QQ*Ka+QQ*Ra)

nxs=nx*4+nxo*3; % size of X
%
nmod=nxs;% number of eigenmodes
MM=zeros(nxs,nxs); % linear matrix
%
% Define dtX
dtKa=zeros(nx,1);
dtRa=zeros(nx,1);
dtA=zeros(nx,1);
dtZ=zeros(nx,1);
dtKo=zeros(nxo,1);
dtRo=zeros(nxo,1);
dtT=zeros(nxo,1);

% Define X
Ka=zeros(nx,1);
Ra=zeros(nx,1);
A=zeros(nx,1);
Z=zeros(nx,1);
Ko=zeros(nxo,1);
Ro=zeros(nxo,1);
T=zeros(nxo,1);% beware is truncated to psi0

% Define Xpert
Xpert=zeros(nxs,1);
perturb=1;

% Index of each field
ixKa=1; ixRa=nx+1; ixA=nx*2+1; ixZ=nx*3+1;% index of variable start
ixKo=nx*4+1; ixRo=nx*4+nxo+1; ixT=nx*4+nxo*2+1; 
nxKa=nx; nxRa=nx*2; nxA=nx*3; nxZ=nx*4;% index of variable end
nxKo=nx*4+nxo; nxRo=nx*4+nxo*2; nxT=nx*4+nxo*3; 

% Start loop on perturbation of X for computation of dtX
for iperturb=1:nxs
%
% Compute Xpertj
Xpert=zeros(nxs,1); Xpert(iperturb,1)=perturb;
Ka=Xpert(ixKa:nxKa);
Ra=Xpert(ixRa:nxRa);
A=Xpert(ixA:nxA);
Z=Xpert(ixZ:nxZ);
Ko=Xpert(ixKo:nxKo);
Ro=Xpert(ixRo:nxRo);
T=Xpert(ixT:nxT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute dtX=MX deterministic model

%%%%%%%%%%%%%
% Compute Eq
Eq=zeros(nx,1); for ix=1:nxo; Eq(ix)=alphaq*T(ix); end 

%%%%%%%%%%%%%
% Compute Eqma with zonal mean removed
Eqzma=Eq-mean(Eq);

%%%%%%%%%%%%%
% Compute dtKa, dtRa 
for ix=2:nx;   
dtKa(ix,1)= - 1/dx*(Ka(ix,1)-Ka(ix-1,1)) -ee*ddu*Ka(ix,1) - HH*A(ix,1)/2; end
dtKa(1,1) = - 1/dx*(Ka(1,1)- Ka(nx,1))   -ee*ddu*Ka(1,1)  - HH*A(1,1)/2;
for ix=1:nx-1; 
dtRa(ix,1) = 1/(3*dx)*(Ra(ix+1,1)-Ra(ix,1)) -ee*ddu*Ra(ix,1) - HH*A(ix,1)/3; end
dtRa(nx,1) = 1/(3*dx)*(Ra(1,1)-Ra(nx,1))    -ee*ddu*Ra(nx,1) - HH*A(nx,1)/3;
      
%%%%%%%%%%%%%
% Compute dtA
for ix=1:nx; dtA(ix,1)=(GG)*MMA(ix,1)*Z(ix,1) + (GG*QQ)*MMA(ix,1)*(Ka(ix,1)+Ra(ix,1)); end

%%%%%%%%%%%%%
% Compute dtZ (Z=Q- QQ*(Ka+Ra))
for ix=1:nx; dtZ(ix,1)=-HH*(1-QQ)*A(ix,1) +muoa*chia*Eqzma(ix,1); end

%%%%%%%%%%%%%
% Compute taux
taux=gammaa*(Ka-Ra);

%%%%%%%%%%%%%
% Compute dtKo,dtRo
for ix=2:nxo;   dtKo(ix,1) = - ee*cc1/dx*(Ko(ix,1)-Ko(ix-1,1))         + ee*muao*chio*cc1/2*taux(ix,1); end
                 dtKo(1,1) = - ee*cc1/dx*(Ko(1,1)-rrw*Ro(1,1))         + ee*muao*chio*cc1/2*taux(1,1);
for ix=1:nxo-1; dtRo(ix,1) =   ee*cc1/(3*dx)*(Ro(ix+1,1)-Ro(ix,1))     - ee*muao*chio*cc1/3*taux(ix,1); end
               dtRo(nxo,1) =   ee*cc1/(3*dx)*(rre*Ko(nxo,1)-Ro(nxo,1)) - ee*muao*chio*cc1/3*taux(nxo,1);

%%%%%%%%%%%%%
% compute H
Hocean=Ko+Ro;

%%%%%%%%%%%%%
% compute dtT
for ix=1:nxo; dtT(ix,1)= -ee*cc1*zeta*Eq(ix,1)+ee*cc1*eta(ix,1)*Hocean(ix,1); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute dtX
dtXX=zeros(nxs,1);
dtXX(ixKa:nxKa,1)=dtKa;
dtXX(ixRa:nxRa,1)=dtRa;
dtXX(ixA:nxA,1)=dtA;
dtXX(ixZ:nxZ,1)=dtZ;
dtXX(ixKo:nxKo,1)=dtKo;
dtXX(ixRo:nxRo,1)=dtRo;
dtXX(ixT:nxT,1)=dtT;
%
% Compute MMij=dtXi/Xpertj
for ixs=1:nxs; MM(ixs,iperturb)=dtXX(ixs,1)/perturb; end
%
end% end of loop on perturbation of X

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear Xpert dtXX
clear dtKa dtRa dtA dtZ dtKo dtRo dtT
clear Ka Ra A Z Ko Ro T
