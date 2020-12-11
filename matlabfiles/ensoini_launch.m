% by Sulian Thual
% 
% main launch program for the ENSO-MJO model
%
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Setup

datafolder='tssgcm'; firstrestart=1; lastrestart=3; % TSS-GCM model
%  datafolder='tssgcmwalker'; firstrestart=1; lastrestart=3; % TSS-GCM model with Walker circulation
%  datafolder='crudeintra'; firstrestart=1; lastrestart=3; % Crude intraseasonal model
%  datafolder='crudeinter'; firstrestart=1; lastrestart=3; % Crude interannual model

ensoini_setup=strcat('ensoini_params_',datafolder,'.m'); % parameter file name
run(ensoini_setup); % run parameter file 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Numerical simulations
%

dolaunch=                  0                 ;% recompute solutions

if dolaunch==1
%
run(ensoini_setup); % run parameter file 

%
for indexrestart=firstrestart:lastrestart% loop over restart files
%
strcat(['restart= ',num2str(indexrestart),' / ',num2str(lastrestart)])% prompt current restart file
run(ensoini_setup); % run parameter file again

ensonum_initial; % compute/read initial conditions for one restart file
ensonum_loop; % update system for one restart file

%
end
%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Analysis

iwind=2; % Starting window for plots

% Model Formulation (profiles etc)
%  ensoplot_formulationprofiles;% profile of sq,so,eta,Gamma dist

% Timeserie and power spectra (for statistics)
%  ensoplot_timemean_efficient;% plots time-mean field (useful for next to remove)
%  ensoplot_timeserie_efficient;% compute/plot timeserie (not too long)
%  ensoplot_powerspectra_efficient;% compute/plot power spectra (very long)

% Hovmoller for an ensemble of restarts(for individual events)
%  ensoplot_hovmoller_efficient;% compute/plot hovmoller very long (omit for next ones)




