% by Sulian Thual
% 
% set initial conditions for each restart of the ENSO-MJO model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case 1: At beginning of a simulation (indexrestart=1), use arbitrary linear solutions as initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if indexrestart==1

% Compute Linear solutions
ensonum_linear; % compute linear matrix
[vects,eigd]=eig(MM); % MM is the linear matrix computed using ensonum_linear.m
eigenvals=zeros(nmod,1);
for ixs=1:nmod; eigenvals(ixs,1)=eigd(ixs,ixs); end
wwp=-imag(eigenvals); % freq
rrp=real(eigenvals); % growth
indexm=1:nmod;[term,indexm]=sort(rrp,'descend');
wwp=wwp(indexm); rrp=rrp(indexm); vects=vects(:,indexm); 
%
% Set Initial Conditions (using linear solutions)
imodshow=1;% index of linear solution used for initiation
iniampl=10; % amplitude of linear solution used for initiation
XXini=iniampl*real(vects(:,imodshow));
clear vects eigd eigenvals wwp rrp indexm

% Set Initial conditions for intraseasonal convection
Abarini=(sq-QQ*so)/(1-QQ)/HH/phi0eq;

% Set time axis
ts=(0:nts-1)*dt*mts; % adim t axis for save
tt=(0:mts*nts-1)*dt; % adim t axis for computing
nt=mts*nts; 



else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case 2: For indexrestart>1, continue simulation (use end of previous restart file as initial conditions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read previous restart file
dfolder=strcat('data_',datafolder); % setup is defined in main
filein=strcat(dfolder,'/ensomjo_',num2str(indexrestart-1), '.nc');
term=ensospe_ncdfgetvar(filein,'Xs'); 

% Set Initial conditions (using previous restart file)
XXini=term(:,end); 

% Set Initial conditions for intraseasonal convection
term=ensospe_ncdfgetvar(filein,'MMAs'); Abarini=term(:,end); 
clear term

% Set time axis for new restart file
term=ensospe_ncdfgetvar(filein,'ts'); tini=term(:,end); 
ts=tini+(1:nts)*dt*mts; % adim t axis (of saves)
tt=(0:mts*nts-1)*dt; % adim t axis for computing
nt=mts*nts; 


end




