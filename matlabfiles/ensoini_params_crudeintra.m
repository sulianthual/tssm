% by Sulian Thual


% all parameters needed to run the ENSO-MJO model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART I: Main Parameters

%%%%%%%%%%%%%%%%%%%%%
% Constants and Dimensionalisation
onehour=3600.;% one hour in seconds
oneday=onehour*24.;% one day in seconds
oneyear=oneday*365; % one year in seconds
onekm=1000.;% one km in meters
ii=complex(0,1);


% Units%%%%%%%%%%%%%
%
% axis
xa=15000.*onekm; %x adim (meters),
ya=1500.*onekm; %y adim (meters)
Ya=330.*onekm; %Y adim (meters)
ta=33*oneday; %tau adim (seconds), !!! (we dont use the t interannual axis in practice)
% atm
oa=1.5; % pot temperature theta adim (K)
qa=1.5; % moisture q adim (K)
ua=5; % velocity adim (m.s-1)
% ocean
Ta=         1.5            ; % dimensional scale SST (K)
Ha=         20.8           ; % dimensional scale h (m)
Ua=         0.25           ; % dimensional scale U ocean (ms-1)
tauxa=      0.009          ;% dimensional scale taux (N.m-2) %XXXXXX%

% Grids%%%%%%%%%%%%%%%%%
%
% x-axis atmosphere
nxfact=1;
nx=64*nxfact; % number of atmospheric points (around equatorial belt)
dx=625*onekm/xa/nxfact;% dx grid step (adim) (instead of 625 for 64 points);
%  nx=65*nxfact; dx=615.3846*onekm/xa/nxfact;% dx grid step (adim) (instead of 625 for 64 points);
xx=(0:nx-1)*dx; % adim x axis
xg=xx*xa/1000/1000; % axis x for graphs, in 1000km
LL=nx*dx;%40000*onekm/xa; % atm size, nondimensional (40 000km)

% x-axis ocean
nxo=28*nxfact;
xxo=(0:nxo-1)*dx; % adim x axis ocean
xgo=xxo*xa/1000/1000; % axis x for graphs, in 1000km
LLo=nxo*dx;%18000*onekm/xa; % ocean size, nondimensional (18 000 km)

% General x-axis
yy=-5:0.1:5; YY=yy;
ny=length(yy);

% y-axis atmosphere (i.e. parabolic cylinder functions)
yyg=yy*ya/1000/1000; % axis y for graphs, in 1000 km
phi0=pi^(-0.25).*exp(-yy.^2/2);% phi0 atm
phi2=(4*pi)^(-0.25)*(2*yy.^2-1).*exp(-yy.^2/2); % phi1 atm
phi0eq=0.7511;% phi0(0) used often in practice

% y-axis ocean
YYg=yy*Ya/1000/1000; % axis Y for graphs, in 1000 km
psi0=pi^(-0.25)*exp(-YY.^2/2);% phi0 atm
psi2=(4*pi)^(-0.25)*(2*YY.^2-1).*exp(-YY.^2/2); % phi1 atm
psi0eq=0.7511;% psi0(0) used often in practice

%t-axis
nts=             3*365*4      ; % length of saved in restart (in days)
mts=             10       ; % number of solved timesteps for one saved timestep
dt=              8*onehour/ta     ; % timestep for solving (adim)(dt=NaN puts dt=dx/2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters (non-dimensional)
% General and Scaling
ee=0.1;% Froude number
cc1=             0.5 ; % c1

% Atmosphere
HH=              22         ;% H in skeleton model
QQ=              0.9        ; % Q in skeleton model
GG=              1.66       ; % Gamma in skeleton model
dda=             10^-8      ; % atm dissipation (unecessary)


% Ocean
eta = 1.5 + 0.5 *tanh(7.5*(xxo-LLo/2)); eta=eta';% eta in ocean 
rrw = 0.5; % boundary reflection west
rre =   1; % boundary reflection east

% SST 
zeta=            8.7        ; % coeff loss from latent heating

% Couplings
alphaq=          0.2042536  ;% latent heat coeff=7*qe*exp(qe*25/1.5)/15, qe=0.0929600
gammaa=           6.529      ; % wind stress coefficient
chia=            0.3086     ; % projection coefficient to atm
chio=            1.3801     ; % projection coefficient to ocean 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters: Additional tunings

simutype=  2; % VERY IMPORTANT, TYPE OF MODEL/SIMULATION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% External sources 
%  sq=          zeros(nx,1)+2.2; % sq in skeleton model at equator
sq=          2.2*(1+0.6*cos(2*pi/LL*xx)'); % warm pool
so=          sq; % so in skeleton model at equator
%  so= 2.2*(1+0.6*cos(2*pi/LL*xx-0.1)');% warm pool walker circulation, phase shift
%  so= 2.2*(1+0.63*cos(2*pi/LL*xx)');% warm pool walker circulation, amplitude shift
MMA=         (sq-QQ*so)/(1-QQ)/HH/phi0eq     ;% mean state A0 (first parabolic cylinder function)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dissipations (unit is interannual, i.e. 33 days time)
ddu=              (1/(30*oneday))*ta*ee; %1e-8      ;% dissipation u-o 
ddu0=             (1/(30*oneday))*ta*ee;% dissipation atm for zonal mean component  
ddq=              (1/(30*oneday))*ta*ee;% dissipation of q
dda=              (1/(30*oneday))*ta*ee;% dissipation of a
ddZ=              (1/(30*oneday))*ta*ee;% dissipation of Z
ddo=              0                 ;%dissipation ocean

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Artificial couplings (ref=1, =0 for no coupling)
muao=             1         ;% coupling coefficient atm to ocean
muoa=             1         ;% coupling coefficient ocean to atm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% External noise sources

sigma_abo=0+zeros(nx,1); % amplitude of perturbation theta
sigma_abZ=0.4+zeros(nx,1); % amplitude of perturbation Z
sigma_aba=0+zeros(nx,1);% amplitude of perturbation a

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gamma distribution
lambda_gsde=1/(30*oneday)*ta*ee; % Gamma distribution relaxation rate
k_sde=2;% shape form of Gamma distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
