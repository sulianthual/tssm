% by Sulian Thual
% 
% plot timeseries from ENSO-MJO model solutions
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot parameters

%

ixshow=1;

dorecomputetime=    0     ;% recompute timeseries
dographtimeserie=  1    ; % do graphs

dodetecttimeserie=   0   ;% FOR TEST: quickly detect interesting events
dotimeserieregular=    1  ;% show regular timeserie
showrestarttimeaxis=   0   ;% show restart file in time axis (DIRTY)
dotimeTHaspectrum =    0;% plot for crude inter (TE, Habar time and spectrum)


showtimeseriepower=   0;% show power spectra
showtimeseriepdf=     1     ;% show pdf
filetimesave=strcat('data_',datafolder,'/timeseries.nc');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read all files

if dorecomputetime==1
% Index of each field
ixKa=1; ixRa=nx+1; ixA=nx*2+1; ixZ=nx*3+1;% index of variable start
ixKo=nx*4+1; ixRo=nx*4+nxo+1; ixT=nx*4+nxo*2+1; 
nxKa=nx; nxRa=nx*2; nxA=nx*3; nxZ=nx*4;% index of variable end
nxKo=nx*4+nxo; nxRo=nx*4+nxo*2; nxT=nx*4+nxo*3; 

% Loop  over restart files
tcon=[0,0];
xconT=[0,0]; % average eastern Pacific
xconTw=[0,0]; % average western Pacific
xconu=[0,0]; % zonal mean Pacific winds intra (average entire Pacific)
xconuw=[0,0]; % zonal mean Pacific winds intra, average western half Pacific
xcona=[0,0]; % one location Ha'
xconma=[0,0];% one location Ham
restartcon=[0,0];% indexofrestart file used

for indexrestart=firstrestart:lastrestart
strcat(['restart= ',num2str(indexrestart),' / ',num2str(lastrestart)])

% Read
dfolder=strcat('data_',datafolder); % setup is defined in main
fileout=strcat(dfolder,'/ensomjo_',num2str(indexrestart), '.nc');

% Concatenate time
term=ensospe_ncdfgetvar(fileout,'ts'); tcon=[tcon,3.3*term];% (dim to days) BEWARE !!!!!
% ts is nondimensional, and time unit is 3.3 days. Spacing of timesteps is 8 hours for saving (0.8 hours for solving). 
pass=length(term); pass=indexrestart+(0:pass-1)/pass; restartcon=[restartcon,pass];%save restart number
%
% Read
XXs=ensospe_ncdfgetvar(fileout,'Xs');
Ka=XXs(ixKa:nxKa,:);
Ra=XXs(ixRa:nxRa,:);
A=XXs(ixA:nxA,:);
T=XXs(ixT:nxT,:);
%  Z=XXs(ixZ:nxZ,:);
%  Ko=XXs(ixKo:nxKo,:);
%  Ro=XXs(ixRo:nxRo,:);
Am=ensospe_ncdfgetvar(fileout,'MMAs');

% Averages
ixmin=15; ixmax=28; % eastern Pacific half
term=Ta*phi0eq*mean(T(ixmin:ixmax,:),1); 
xconT=[xconT,term];  % concatenate

ixmin=1; ixmax=13; % western Pacific half
term=Ta*phi0eq*mean(T(ixmin:ixmax,:),1); 
xconTw=[xconTw,term];  % concatenate


ixmin=1; ixmax=28; % average entire Pacific
pass=Ka(ixmin:ixmax,:)-Ra(ixmin:ixmax,:);
term=ua*phi0eq*mean(pass,1); % one point
xconu=[xconu,term];  % concatenate

ixmin=1; ixmax=14; % average entire Pacific
pass=Ka(ixmin:ixmax,:)-Ra(ixmin:ixmax,:);
term=ua*phi0eq*mean(pass,1); % one point
xconuw=[xconuw,term];  % concatenate

term=HH*oa/(ta/oneday)*phi0eq*squeeze(A(ixshow,:)); % one point
xcona=[xcona,term];  % concatenate

term=HH*oa/(ta/oneday)*phi0eq*squeeze(Am(ixshow,:)); % one point
xconma=[xconma,term];  % concatenate


end

% Take 
tcon=tcon(3:end);
restartcon=restartcon(3:end); restartcon=round(restartcon);
xconT=xconT(3:end); 
xconTw=xconTw(3:end); 
xcona=xcona(3:end); 
xconu=xconu(3:end); 
xconuw=xconuw(3:end); 
xconma=xconma(3:end); 

% save
ensospe_ncdfmakevar(filetimesave,'Te',{'one','tt'},xconT,NaN,2); 
ensospe_ncdfmakevar(filetimesave,'Tw',{'one','tt'},xconTw,NaN,1); 
ensospe_ncdfmakevar(filetimesave,'uzm',{'one','tt'},xconu,NaN,1);
ensospe_ncdfmakevar(filetimesave,'uw',{'one','tt'},xconuw,NaN,1);
ensospe_ncdfmakevar(filetimesave,'Ha',{'one','tt'},xcona,NaN,1);
ensospe_ncdfmakevar(filetimesave,'Ham',{'one','tt'},xconma,NaN,1);
ensospe_ncdfmakevar(filetimesave,'ts',{'one','tt'},tcon,NaN,1); 
ensospe_ncdfmakevar(filetimesave,'irestart',{'one','tt'},restartcon,NaN,1);

end% end redocomputation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All graphs
if dographtimeserie==1
% Read
xconT=ensospe_ncdfgetvar(filetimesave,'Te');
xconTw=ensospe_ncdfgetvar(filetimesave,'Tw');
xconu=ensospe_ncdfgetvar(filetimesave,'uzm');
xconuw=ensospe_ncdfgetvar(filetimesave,'uw');
xcona=ensospe_ncdfgetvar(filetimesave,'Ha');
xconma=ensospe_ncdfgetvar(filetimesave,'Ham');
tcon=ensospe_ncdfgetvar(filetimesave,'ts');% is now in days cf correction from raw outputs !!!!(was adim)
restartcon=ensospe_ncdfgetvar(filetimesave,'irestart'); 
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Quick detection method 

if dodetecttimeserie==1
figure(iwind); iwind=iwind+1; clf;
ix=1; iy=5; icount=1;

ndcount=100;%100;% 1/5 of total simulation time
for icountc=0:4;% divide in 5 plots
subplot(iy,ix,icountc+1);
ixminc=icountc*ndcount; ixmaxc=icountc*ndcount+ndcount;
if icountc==0; ixminc=2; end
if icountc==4; ixmaxc=5*ndcount; end
ixminc=ensospe_searchclosest(restartcon,ixminc); ixmaxc=ensospe_searchclosest(restartcon,ixmaxc);

plot(restartcon(ixminc:ixmaxc),xconT(ixminc:ixmaxc),'k-')
hold on; plot(restartcon(ixminc:ixmaxc),restartcon(ixminc:ixmaxc)*0+5,'r-');
ylim([-10 10]);
if icountc==0; title('TE detection'); end
if icountc==4; xlabel('restart file'); end
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regular Plots
if dotimeserieregular==1
figure(iwind); iwind=iwind+1; clf;

if showrestarttimeaxis==1; tcon=restartcon*365; end;% DIRTY QUICK FIX

ix=1; iy=4; icount=1;

xrange=[min(tcon/365),max(tcon/365)];
%  xrange=[510 530];
%  xrange=[0 2];
addlinefilt=1;

if 1==1
subplot(iy,ix,1); icount=icount+1;
plot(tcon/365,xconT)
if addlinefilt==1
hold on; plot(tcon/365,movmean(xconT,365),'r-')
%  hold on; plot(tcon/365,Ta*movmean(xconT,[30 0]),'k-');% backward 30 days
%  hold on; plot(tcon/365,Ta*movmean(xconT,[90 0]),'r-');% backward 90 days
end
%  xlabel('time (years)');
xlim(xrange)
title('$T_{E}(K)$','interpreter','latex');
%  ylim([-7 2])
%  ylim([-10 1])
end

if 1==1
subplot(iy,ix,2); icount=icount+1;
plot(tcon/365,xconu)
if addlinefilt==1
hold on; plot(tcon/365,movmean(xconu,365),'r-')
%  hold on; plot(tcon/365,ua*movmean(xconu,[30 0]),'k-');% backward 30 days
%  hold on; plot(tcon/365,ua*movmean(xconu,[90 0]),'r-');% backward 90 days
end
xlabel('time (years)','interpreter','latex');
xlim(xrange)
title('$u_{zm}(m.s^{-1})$','interpreter','latex');
ylim([-25 25])
end

if 1==1
subplot(iy,ix,3); icount=icount+1;
plot(tcon/365,xcona)
if addlinefilt==1
hold on; plot(tcon/365,movmean(xcona,365),'r-')
%  hold on; plot(tcon/365,HH*oa/(ta/oneday)*movmean(xcona,[30 0]),'k-')% backward 30 days
%  hold on; plot(tcon/365,HH*oa/(ta/oneday)*movmean(xcona,[90 0]),'r-')% backward 90 days
end
%  xlabel('time (years)');
xlim(xrange)
title('$\overline{H}a (K.day^{-1})$','interpreter','latex');
%  ylim([-0.1 7])
ylim([-0.1 5])
%   ylim([-3 3])
end

if 1==1
subplot(iy,ix,4); icount=icount+1;
plot(tcon/365,xconma)
if addlinefilt==1
hold on; plot(tcon/365,movmean(xconma,365),'r-')
%  hold on; plot(tcon/365,HH*oa/(ta/oneday)*movmean(xconma,[30 0]),'k-')% backward 30 days
%  hold on; plot(tcon/365,HH*oa/(ta/oneday)*movmean(xconma,[90 0]),'r-')% backward 90 days
end
xlabel('time (years)','interpreter','latex');
xlim(xrange)
title('$\overline{H}\overline{a} (K.day^{-1})$','interpreter','latex');
%  ylim([0.1 0.4])
end

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% power spectra TE
if showtimeseriepower==1
varplot=xconT; % variable to plot
figure(iwind); iwind=iwind+1;, clf;
T=tcon(end)/365-tcon(1)/365; % time of interval (years): REVOIR
N = length(tcon);% number of points
p=abs(fft(varplot)); % absolute value of fft
p=p(1:N/2).^2; % take power of positive freq half
p=p./N^2; % normalise (note we should normalise by (N/2)^2 if cutting Nyquist)
freq=[0:N/2-1]/T; % find corresponding frequency in cpd (missing N/2 ?)
%plot(freq,p); % power
semilogy(freq,p); % plot on semilog scale
hold on; plot([1,1],[10^-12,10^3],'k'); 
%  hold on; plot([0.5,0.5],[10^-12,10^3],'r'); 
%  hold on; plot([0.33,0.33],[10^-12,10^3],'r'); 
%  hold on; plot([0.25,0.25],[10^-12,10^3],'k'); 
%  hold on; plot([0.2,0.2],[10^-12,10^3],'r'); 
xlim([0,2])
ylim([10^-8,10^1])
xlabel('$Frequency (yr^{-1})$','interpreter','latex');
title('Power Spectra $T_{E}$','interpreter','latex');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % pdf
 
 
 if showtimeseriepdf==1
 
 
figure(31); clf; 
%%%%
xconsh=xconT;
[fi,xx]=ksdensity(xconsh);%fi = ensospe_smoothn(fi,[1,10]);

mu = mean(xconsh);
sig = var(xconsh);
skw=skewness(xconsh);
fi_g = normpdf(xx,mu,sqrt(sig));
hold on
plot(xx,fi,'linewidth',2)
plot(xx,fi_g,'--r','linewidth',2)
set(gca,'fontsize',12)
title('$T_{E}$','interpreter','latex')
%  xlim(xrange)
%  ylabel(datafolder);
xlabel({strcat('m=',num2str(mu));strcat('v=',num2str(sig));strcat('s=',num2str(skw))});
ylim([0 max(fi)*1.1])
box on
 
 
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
% special graph crude intraseasonal : timeseries of TE, HAbarw, spectrum and pdf
 if dotimeTHaspectrum==1

figure(iwind); iwind=iwind+1; clf;

ix=2; iy=4; icount=1;

%  xrange=[min(tcon/365),max(tcon/365)];
xrange=[0 200]+1600;
%  xrange=[0 2];
addlinefilt=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT TE
if 1==1
subplot(iy,ix,[1,2]); icount=icount+1;
plot(tcon/365,xconT)
hold on; plot(tcon/365,tcon*0,'k-')
if addlinefilt==1
%  hold on; plot(tcon/365,movmean(xconT,365),'r-')
%  hold on; plot(tcon/365,Ta*movmean(xconT,[30 0]),'k-');% backward 30 days
%  hold on; plot(tcon/365,Ta*movmean(xconT,[90 0]),'r-');% backward 90 days
end
%  xlabel('time (years)','interpreter','latex');
xlim(xrange)
title('(a) $T_{E}$','interpreter','latex');
ylabel('$K$','interpreter','latex');
ylim([-7 7])
%  ylim([-10 1])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT HABAR
if 1==1
subplot(iy,ix,[3,4]); icount=icount+1;
plot(tcon/365,xconma)
if addlinefilt==1
hold on; plot(tcon/365,movmean(xconma,365),'r-')
%  hold on; plot(tcon/365,HH*oa/(ta/oneday)*movmean(xconma,[30 0]),'k-')% backward 30 days
%  hold on; plot(tcon/365,HH*oa/(ta/oneday)*movmean(xconma,[90 0]),'r-')% backward 90 days
end
xlabel('time (years)','interpreter','latex');
xlim(xrange)
title('(b) $\overline{H}\overline{a}$','interpreter','latex');
ylabel('$K.day^{-1}$','interpreter','latex');
ylim([0 1])
end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
 % power spectrum TE
 
subplot(iy,ix,5); icount=icount+1;
 
dologspectrumT=1;
 
varplot=xconT; % variable to plot
T=(tcon(end)-tcon(1))/365; % time of interval (days)
N = length(tcon);% number of points
pow=abs(fft(varplot)); % absolute value of fft
pow=pow(1:N/2).^2; % take power of positive freq half
pow=pow./N^2; % normalise (note we should normalise by (N/2)^2 if cutting Nyquist)
freq=[0:N/2-1]/T; % find corresponding frequency in cpd (missing N/2 ?)
%  pow=log10(pow);


% Plot full
if dologspectrumT==1
semilogy(pow,freq);
loglog(pow,freq);
else
plot(pow,freq);
end
%  semilogx(pow,freq);

% Compute extremely smoothed version
ntsmoo=     101; % t size of smooth (odd) (ref=25)
nsmoo=5;% number of times smoothing

for ismoo=1:nsmoo
pow=exp(ensospe_smoothn(log(pow),[1,ntsmoo]));% smooth
end
hold on; plot(pow,freq,'r-','linewidth',2);

%  hold on; plot(pow,freq,'r-','linewidth',2);
%
yticks([0.01 0.1 1 10]);
xrange=[10e-15 1];
%  xticks([10e-16 10e-11 1]);
hold on; plot(xrange,xrange*0+1/30*365,'k--','linewidth',2);
hold on; plot(xrange,xrange*0+1/90*365,'k--','linewidth',2);
%  hold on; plot(xrange,xrange*0+1/10*365,'k-','linewidth',2);
%  hold on; plot(xrange,xrange*0+1/10,'k-','linewidth',2);
xlim(xrange)
if dologspectrumT==1
ylim([0.08,365/8]); % intra; goes from 10 years to 5 days
else
ylim([0,15]); % intra
end
%  ylim([0,0.01]); % inter
%
ylabel('$yr^{-1}$','interpreter','latex');
xlabel('$K^{2}$','interpreter','latex');
%  ylabel('$Frequency (cpd)$','interpreter','latex');
title('(c) Spectrum: $T_{E}$','interpreter','latex');
%  colorbar('location','southoutside','visible','off');
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
 % power spectrum Habar
 
subplot(iy,ix,6); icount=icount+1;
 
dologspectrumT=1;
 
varplot=xconma-mean(xconma(:)); % variable to plot
T=(tcon(end)-tcon(1))/365; % time of interval (days)
N = length(tcon);% number of points
pow=abs(fft(varplot)); % absolute value of fft
pow=pow(1:N/2).^2; % take power of positive freq half
pow=pow./N^2; % normalise (note we should normalise by (N/2)^2 if cutting Nyquist)
freq=[0:N/2-1]/T; % find corresponding frequency in cpd (missing N/2 ?)
%  pow=log10(pow);


% Plot full
if dologspectrumT==1
semilogy(pow,freq);
loglog(pow,freq);
else
plot(pow,freq);
end
%  semilogx(pow,freq);

% Compute extremely smoothed version
ntsmoo=     101; % t size of smooth (odd) (ref=25)
nsmoo=5;% number of times smoothing

for ismoo=1:nsmoo
pow=exp(ensospe_smoothn(log(pow),[1,ntsmoo]));% smooth
end
hold on; plot(pow,freq,'r-','linewidth',2);

%  hold on; plot(pow,freq,'r-','linewidth',2);
%
yticks([0.01 0.1 1 10]);
xrange=[10e-15 1];
%  xticks([10e-16 10e-11 1]);
hold on; plot(xrange,xrange*0+1/30*365,'k--','linewidth',2);
hold on; plot(xrange,xrange*0+1/90*365,'k--','linewidth',2);
%  hold on; plot(xrange,xrange*0+1/10*365,'k-','linewidth',2);
%  hold on; plot(xrange,xrange*0+1/10,'k-','linewidth',2);
xlim(xrange)
if dologspectrumT==1
ylim([0.08,365/8]); % intra; goes from 10 years to 5 days
else
ylim([0,15]); % intra
end
%  ylim([0,0.01]); % inter
%
ylabel('$yr^{-1}$','interpreter','latex');
xlabel('$K^{2}.day^{-2}$','interpreter','latex');
%  ylabel('$Frequency (cpd)$','interpreter','latex');
title('(d) Spectrum: $\overline{H} \overline{a}$','interpreter','latex');
%  colorbar('location','southoutside','visible','off');
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Pdf TE
 
 subplot(iy,ix,7);
 xconsh=xconT;
[fi,xx]=ksdensity(xconsh);%fi = ensospe_smoothn(fi,[1,10]);

mu = mean(xconsh);
sig = var(xconsh);
skw=skewness(xconsh);
fi_g = normpdf(xx,mu,sqrt(sig));
hold on
plot(xx,fi,'linewidth',2)
plot(xx,fi_g,'--r','linewidth',2)
%  set(gca,'fontsize',12)
title('(e) Pdf:$T_{E}$','interpreter','latex')
xlim([-15 15])
%  ylabel(datafolder);
%  xlabel({strcat('m=',num2str(mu));strcat('v=',num2str(sig));strcat('s=',num2str(skw))});
xlabel('$K$','interpreter','latex')
ylim([0 max(fi)*1.1])
box on

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Pdf Habar
 subplot(iy,ix,8);
 xconsh=xconma;
[fi,xx]=ksdensity(xconsh);%fi = ensospe_smoothn(fi,[1,10]);

mu = mean(xconsh);
sig = var(xconsh);
skw=skewness(xconsh);
%  fi_g = normpdf(xx,mu,sqrt(sig));
hold on
plot(xx,fi,'linewidth',2)
%  plot(xx,fi_g,'--r','linewidth',2)
%  set(gca,'fontsize',12)
title('(f) Pdf: $\overline{H} \overline{a}$','interpreter','latex')
xlabel('$K.day^{-1}$','interpreter','latex')
xlim([0 0.6])
%  ylabel(datafolder);
%  xlabel({strcat('m=',num2str(round(mu,1)));strcat('v=',num2str(round(sig,2)));strcat('s=',num2str(round(skw,2)))});
ylim([0 max(fi)*1.1])
box on


 
 end% plot crude inter
%%%%%%%%%%%%%%%%%%%%%%%%%%
end% dographtimeserie
%%%%%%%%%%%%%%%%%%%%%%%%%
