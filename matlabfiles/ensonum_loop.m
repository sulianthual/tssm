% by Sulian Thual
% 
% run the ENSO-MJO model over one restart
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving
%
%
if 1==1
%
% Output file
dfolder=strcat('data_',datafolder); % setup is defined in main
fileout=strcat(dfolder,'/ensomjo_',num2str(indexrestart), '.nc');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation of Run
%
% Prepare Netcdf (Saved Variables): we store the entire state X
term=zeros(nxs,nts); term(:,1)=XXini;
ensospe_ncdfmakevar(fileout,'Xs',{'X','T'},term,NaN,2);
ensospe_ncdfmakevar(fileout,'ts',{'one','T'},ts,NaN,1); term=0;% BEWARE, ts is nondim ! save is each dt*mts=3.3 days !!!!!
term=zeros(nx,nts); term(:,1)=Abarini;
ensospe_ncdfmakevar(fileout,'MMAs',{'X1','T'},term,NaN,1);% add A_bar
%Initialise run and increment variables
XX=XXini; XXe=XX*0; Abar=Abarini; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time loop of Run
%
ixKa=1; ixRa=nx+1; ixA=nx*2+1; ixZ=nx*3+1;% index of variable start
ixKo=nx*4+1; ixRo=nx*4+nxo+1; ixT=nx*4+nxo*2+1; 
nxKa=nx; nxRa=nx*2; nxA=nx*3; nxZ=nx*4;% index of variable end
nxKo=nx*4+nxo; nxRo=nx*4+nxo*2; nxT=nx*4+nxo*3; 

%
for its=1:nts
for it=1:mts

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Linear component 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART I: Solve background mean flow: Ko, Ro, T (with MMA)

%%%%%%%%%%%
% Part IA: Compute background convective activity Abar
%          IMPORTANT: impose >=0 if Eqzma too negative

Eq=zeros(nx,1); Eq(1:nxo)=alphaq*XXe(ixT:nxT,1); 
Eqzma=Eq-mean(Eq);

if simutype==1; % perturbed Abar, Gamma distribution
Abar0=1/(HH*(1-QQ))*(chia*Eqzma+sq/phi0eq-QQ*so/phi0eq);
Abar0(Abar0<=0)=10e-5;
dtAbar1= -lambda_gsde*(Abar-Abar0); % relaxation term
dtAbar2= sqrt(2/k_sde*lambda_gsde*Abar.*Abar0).*randn(nx,1);
Abar=Abar+dt*dtAbar1 +sqrt(dt)*dtAbar2; 

else

Abar=1/(HH*(1-QQ))*(chia*Eqzma+sq/phi0eq-QQ*so/phi0eq);
Abar(Abar<=0)=10e-5;% impose Abar>0

end
clear Eq Eqzma


%%%%%%%%%%%%%
% Part IB: compute interannual atm
term=so/phi0eq-HH*Abar;
tauxinter=gammaa*ensonum_oldatm(3,nx,dx,10e-8,1, (3/2)*term); 
clear term
%  tauxinter=tauxinter*0; % TEST: put interannual atm to zero

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part II: Solve intraseasonal atmosphere (all anomalies)


%%%%%%%%%%%%%
% Part IIA: Solve KA, RA (optional additive noise)
%
Abo=1/sqrt(dt)*sigma_abo.*randn(nx,1); % perturbation
term=-1/2*( HH*XX(ixA:nxA,1) +Abo );
XXe(ixKa:nxKa,1)=ensonum_wavefourier(XX(ixKa:nxKa,1),ddu,ddu,1,term,nx,dx,dt);
term=-1/3*( HH*XX(ixA:nxA,1) +Abo );
XXe(ixRa:nxRa,1)=ensonum_wavefourier(XX(ixRa:nxRa,1),ddu,ddu,-1/3,term,nx,dx,dt);
clear term Abo 


%%%%%%%%%%%%%
% Part IIB: Solve Z (optional additive noise)
%
AbZ=1/sqrt(dt)*sigma_abZ.*randn(nx,1); % perturbation
dtZ=-HH*(1-QQ)*( XX(ixA:nxA,1)+AbZ ) -ddZ*XX(ixZ:nxZ,1);
XXe(ixZ:nxZ,1)=XX(ixZ:nxZ,1)+dt*dtZ;
clear dtZ AbZ

%%%%%%%%%%%%
% Part IIIB: Solve A (linear/nonlinear-gamma, optional noise)

% Part IIIB1: Solve A deterministic part

if simutype==2; % REFERENCE: A linearized (optional additive noise)
%
Aba=1/sqrt(dt)*sigma_aba.*randn(nx,1); % perturbation
dtA=GG*Abar.*XX(ixZ:nxZ,1) + (GG*QQ)*Abar.*( XXe(ixKa:nxKa,1)+XXe(ixRa:nxRa,1) ) -dda*XX(ixA:nxA,1);% linearized 
XXe(ixA:nxA,1)=XX(ixA:nxA,1)+dt*dtA+dt*Aba;
clear dtA Aba
end


if simutype==3; % REFERENCE: A nonlinear, gamma distribution
%
Atot=XX(ixA:nxA,1)+ Abar; 
Atot(Atot<=0)=10e-5;% ensure total A positive
dtA=GG*Atot.*XX(ixZ:nxZ,1) + GG*QQ*Atot.*( XXe(ixKa:nxKa,1)+XXe(ixRa:nxRa,1) );
dtAdiss=-lambda_gsde*XX(ixA:nxA,1);
dtAgamma=sqrt(2/k_sde*lambda_gsde*Atot.*Abar).*randn(nx,1)    ; % 
%  dtAgamma=sqrt(lambda_gsde*Atot.*Abar).*randn(nx,1)    ; % old ERROR, missing ee
%
Ae=XX(ixA:nxA,1)+dt*dtA+dt*dtAdiss+sqrt(dt)*dtAgamma;
Ae(Ae<=-Abar)=-Abar(Ae<=-Abar)+10e-5;% ensure total A positive
XXe(ixA:nxA,1)=Ae; 
%
clear dtA dtAdiss dtAgamma Ae
end

%%%%%%%%%%%%
if simutype==1;% TEST, put all deterministic intraseasonal atm to rest
XXe(ixKa:nxKa,1)=0; 
XXe(ixRa:nxRa,1)=0; 
XXe(ixZ:nxZ,1)=0; 
XXe(ixA:nxA,1)=0; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part III: Solve ocean

%%%%%%%%%%%%%
% Part IIIA: Compute dtKo,dtRo
dtKo=zeros(nxo,1); dtRo=zeros(nxo,1);
taux=tauxinter+gammaa*(XXe(ixKa:nxKa)-XXe(ixRa:nxRa)); 

Ko=XX(ixKo:nxKo); Ro=XX(ixRo:nxRo); 
for ix=2:nxo;   dtKo(ix,1) = - ee*cc1/dx*(Ko(ix,1)-Ko(ix-1,1))         + ee*muao*chio*cc1/2*taux(ix,1); end
                 dtKo(1,1) = - ee*cc1/dx*(Ko(1,1)-rrw*Ro(1,1))         + ee*muao*chio*cc1/2*taux(1,1);
for ix=1:nxo-1; dtRo(ix,1) =   ee*cc1/(3*dx)*(Ro(ix+1,1)-Ro(ix,1))     - ee*muao*chio*cc1/3*taux(ix,1); end
               dtRo(nxo,1) =   ee*cc1/(3*dx)*(rre*Ko(nxo,1)-Ro(nxo,1)) - ee*muao*chio*cc1/3*taux(nxo,1);
XXe(ixKo:nxKo,1)=Ko+dt*dtKo;
XXe(ixRo:nxRo,1)=Ro+dt*dtRo;
clear Ko Ro dtKo dtRo tauxinter

%%%%%%%%%%%%%
% Part IIIB: compute dtT
dtT= -ee*cc1*zeta*alphaq*XX(ixT:nxT,1)+ee*cc1*eta.*( XXe(ixKo:nxKo,1)+XXe(ixRo:nxRo,1) ); 
XXe(ixT:nxT,1)=XX(ixT:nxT,1)+dt*dtT;
clear dtT

if 0==1;% TEST: put ocean to rest
XXe(ixKo:nxKo,1)=XXe(ixKo:nxKo,1)*0;
XXe(ixRo:nxRo,1)=XXe(ixRo:nxRo,1)*0;
XXe(ixT:nxT,1)=XXe(ixT:nxT,1)*0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%% Increment variables
XX=XXe;

end%loop it

%%%%% Save 
ensospe_ncdfmakevar(fileout,'Xs',{'X','T'},XXe,[1,its],0);
ensospe_ncdfmakevar(fileout,'MMAs',{'X','T'},Abar,[1,its],0);
end;% loop its
%
end; % 1==1 key for run
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
