% by Sulian Thual
% 
% contour power spectra in the atmosphere
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%
dorecomputepower=  0  ;% RECOMPUTE POWER (VERY LONG !!!!)
dographpower=    1   ;% do graph power kw
dographspectrum= 1   ;% do graph Te, uw spectrum
%

%
filepowersave=strcat('data_',datafolder,'/powerspectra.nc');
if dorecomputepower==1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read all files

% Index of each field
ixKa=1; ixRa=nx+1; ixA=nx*2+1; ixZ=nx*3+1;% index of variable start
ixKo=nx*4+1; ixRo=nx*4+nxo+1; ixT=nx*4+nxo*2+1; 
nxKa=nx; nxRa=nx*2; nxA=nx*3; nxZ=nx*4;% index of variable end
nxKo=nx*4+nxo; nxRo=nx*4+nxo*2; nxT=nx*4+nxo*3; 

% Loop  over variable/restart restart files
for ivar=1:4
xcon=((1:nx)*0)';% term to concatenate (redefine for each variable)
for indexrestart=firstrestart+1:lastrestart
strcat(['var= ', num2str(ivar),', restart= ',num2str(indexrestart),' / ',num2str(lastrestart)])

% Read
dfolder=strcat('data_',datafolder); % setup is defined in main
fileout=strcat(dfolder,'/ensomjo_',num2str(indexrestart), '.nc');

% Concatenate variables
XXs=ensospe_ncdfgetvar(fileout,'Xs');

if ivar==1; 
Ka=XXs(ixKa:nxKa,:);
Ra=XXs(ixRa:nxRa,:);
xread=Ka-Ra; 
end
if ivar==2; 
Ka=XXs(ixKa:nxKa,:);
Ra=XXs(ixRa:nxRa,:);
xread=-Ka-Ra;  
end
if ivar==3; 
Ka=XXs(ixKa:nxKa,:);
Ra=XXs(ixRa:nxRa,:);
Z=XXs(ixZ:nxZ,:);
xread=Z+QQ*(Ka+Ra); 
end
if ivar==4; 
xread=XXs(ixA:nxA,:);
end

xcon=[xcon,xread];  % concatenate

end% loop restarts

% Remove first dummy
xcon=xcon(:,2:end); 


% remove temporal mean
for ix=1:nx
xcon(ix,:)=xcon(ix,:)-mean(xcon(ix,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute power

term=xcon; 
clear xcon;
pass=size(term); nth=pass(2);
term=fftshift(fft2(term)); 
term=abs(term).^2/(nx*nth)^2; 
kg=ensospe_fftkspe(nx,dx); % adim
wg=ensospe_fftkspe(nth,dt*mts); %adim
kg=kg/xa*40000*1000/(2*pi);
wg=wg/ta*oneday; 
kg=-kg; % why ? this is no longer due to transpositions '(instead of .') of complex matrix ?
term=log10(term);
%
% Further treatment: restrict w range
wgmin=ensospe_searchclosest(wg,-0.1); wgmax=ensospe_searchclosest(wg,0.2);
term=term(:,wgmin:wgmax); 
wg=wg(wgmin:wgmax);


% Further treatment: smooth wg
nwsmoo=     101; % wg size of smooth (odd) (ref=100)
term=ensospe_smoothn(term,[1,nwsmoo]);% smooth

% Further treatment: undersample
nwu=101; % w-space between samples 
term=term(:,1:nwu:end); 
wg=wg(1:nwu:end);

% Write down
if ivar==1; ensospe_ncdfmakevar(filepowersave,'u',{'k','w'},term,NaN,2); 
ensospe_ncdfmakevar(filepowersave,'kg',{'one','k'},kg,NaN,1)
ensospe_ncdfmakevar(filepowersave,'wg',{'one','w'},wg,NaN,1)
end
if ivar==2; ensospe_ncdfmakevar(filepowersave,'o',{'k','w'},term,NaN,1); end
if ivar==3; ensospe_ncdfmakevar(filepowersave,'q',{'k','w'},term,NaN,1); end
if ivar==4; ensospe_ncdfmakevar(filepowersave,'a',{'k','w'},term,NaN,1); end
clear term

end% loop ivar
%  stop;% just stop here
end% REDO COMPUTATION




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part analysis

if dographpower==1

%  colormap('default');
colormap('jet');
dolevs=1;
if 1==1; % SETUP FOR MAIN WP SIMULATION
levelsu=[0:0.05:1]*3.3-9; 
levelso=[0:0.05:1]*7.1-12; 
levelsq=[0:0.05:1]*3-9; 
levelsa=[0:0.05:1]*4.6-12; 
end

if 1==1;% SETUP FOR CRUDE INTRA
levelsu=[0:0.05:1]*6.5-11; 
levelso=[0:0.05:1]*8.5-13; 
levelsq=[0:0.05:1]*3-9; 
levelsa=[0:0.05:1]*8-14; 
end

if 0==1; % SETUP FOR MAIN WP SIMULATION
levelsu=[0:0.05:1]*4.3-9; 
levelso=[0:0.05:1]*8.5-13; 
levelsq=[0:0.05:1]*3-9; 
levelsa=[0:0.05:1]*8-14; 
end

xconu=ensospe_ncdfgetvar(filepowersave,'u');
xcono=ensospe_ncdfgetvar(filepowersave,'o');
xconQ=ensospe_ncdfgetvar(filepowersave,'q');
xconA=ensospe_ncdfgetvar(filepowersave,'a');
kg=ensospe_ncdfgetvar(filepowersave,'kg');
wg=ensospe_ncdfgetvar(filepowersave,'wg');

% Further treatment: restrict k,w range
[pass,iksort]=sort(kg); % reput as increasing
xconu=xconu(iksort,:);
xcono=xcono(iksort,:);
xconQ=xconQ(iksort,:);
xconA=xconA(iksort,:);
kg=kg(iksort);

wgmin=ensospe_searchclosest(wg,-0.01); wgmax=ensospe_searchclosest(wg,0.1);
kgmin=ensospe_searchclosest(kg,-5); kgmax=ensospe_searchclosest(kg,5);
xconu=xconu(kgmin:kgmax,wgmin:wgmax); 
xcono=xcono(kgmin:kgmax,wgmin:wgmax); 
xconQ=xconQ(kgmin:kgmax,wgmin:wgmax); 
xconA=xconA(kgmin:kgmax,wgmin:wgmax); 
kg=kg(kgmin:kgmax);
wg=wg(wgmin:wgmax);

% Further treatment: smooth
nxsmoo=     1; % x size of smooth (odd) (ref=1)
ntsmoo=     11; % t size of smooth (odd) (ref=25)
nsmoo=100;% number of times smoothing
%  nxsmoo=1; ntsmoo=1; nsmoo=1;
nxsmoo=1; ntsmoo=11; nsmoo=5;
for ismoo=1:nsmoo
xconu=ensospe_smoothn(xconu,[nxsmoo,ntsmoo]);% smooth
xcono=ensospe_smoothn(xcono,[nxsmoo,ntsmoo]);% smooth
xconQ=ensospe_smoothn(xconQ,[nxsmoo,ntsmoo]);% smooth
xconA=ensospe_smoothn(xconA,[nxsmoo,ntsmoo]);% smooth
end
%

%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contour

% restrict range

%  xxrange=[-10,10]; yyrange=[0,0.2];

addlineintra=1; % add 30 and 90 days
adddisp=1; % add dots dispersion relation

%%%%%%%%%%%%%%%%%

% add dispersion curves linear stability
if adddisp==1
% Parameters
ekr=0.1*(1/(30*oneday))*ta; % assumed dissip all variables, must have e=0.1 cf dim time=33 days
sorefdisp=2.2; % assumed sq/so here

nkk=length(kg);
wgrdisp=zeros(nkk,4);
wgidisp=zeros(nkk,4);

for ikk=1:length(kg)
kkhere=kg(ikk)*xa/40000/1000*(2*pi); % re-nondim
[wr,wi,vects]=ensospe_stabskelnew(kkhere,ekr,QQ,GG,HH,sorefdisp);
wgrdisp(ikk,:)=wr;
wgidisp(ikk,:)=wi;
end
wgrdisp=wgrdisp/ta*oneday;

end



%%%%%%%%%%%%%%%%%
figure(iwind); clf; iwind=iwind+1;
ix=3; iy=2; icount=1;
xxrange=[-5,5]; yyrange=[0,0.1];
colormap(jet)



%%%%%%%%%%%%%%%%%%%
% PART 1: POWER SPECTRA kw

% u
subplot(iy,ix,1); 
term=xconu; levels=levelsu;
if dolevs==1
contourf(kg,wg,term',levels, 'Linestyle', 'none'); 
caxis([levels(1) levels(end)]);% VERY IMPORTANT !!!!
colorbar('southoutside')
else
contourf(kg,wg,term','LineStyle','none');
end
xlim(xxrange); ylim(yyrange)
colorbar('southoutside')
xlabel('wavenumber(2$\pi$/40,000km)','interpreter','latex')
ylabel('frequency(cpd)','interpreter','latex')
title('(a) $u^{\prime}$','interpreter','latex');
if addlineintra==1; hold on; plot(wg*0,wg,'k-','linewidth',2); end
if addlineintra==1; hold on; plot(kg,kg*0+1/30,'k--','linewidth',2); end
if addlineintra==1; hold on; plot(kg,kg*0+1/90,'k--','linewidth',2); end
if adddisp==1; for idsp=1:4; hold on; plot(kg,wgrdisp(:,idsp),'ko'); end; end

% A
subplot(iy,ix,2); 
term=xconA; levels=levelsa;
if dolevs==1
contourf(kg,wg,term',levels, 'Linestyle', 'none'); 
caxis([levels(1) levels(end)]);% VERY IMPORTANT !!!!
colorbar('southoutside')
else
contourf(kg,wg,term','LineStyle','none');
end
xlim(xxrange); ylim(yyrange)
colorbar('southoutside')
xlabel('wavenumber(2$\pi$/40,000km)','interpreter','latex')
%  ylabel('frequency(cpd)','interpreter','latex')
title('(b) $a^{\prime}$','interpreter','latex');
if addlineintra==1; hold on; plot(wg*0,wg,'k-','linewidth',2); end
if addlineintra==1; hold on; plot(kg,kg*0+1/30,'k--','linewidth',2); end
if addlineintra==1; hold on; plot(kg,kg*0+1/90,'k--','linewidth',2); end
if adddisp==1; for idsp=1:4; hold on; plot(kg,wgrdisp(:,idsp),'ko'); end; end

% theta
subplot(iy,ix,4); icount=icount+1;
term=xcono; levels=levelso;
if dolevs==1
contourf(kg,wg,term',levels, 'Linestyle', 'none'); 
caxis([levels(1) levels(end)]);% VERY IMPORTANT !!!!
colorbar('southoutside')
else
contourf(kg,wg,term','LineStyle','none');
end
xlim(xxrange); ylim(yyrange)
colorbar('southoutside')
xlabel('wavenumber(2$\pi$/40,000km)','interpreter','latex')
ylabel('frequency(cpd)','interpreter','latex')
title('(d) $\theta^{\prime}$','interpreter','latex');
if addlineintra==1; hold on; plot(wg*0,wg,'k-','linewidth',2); end
if addlineintra==1; hold on; plot(kg,kg*0+1/30,'k--','linewidth',2); end
if addlineintra==1; hold on; plot(kg,kg*0+1/90,'k--','linewidth',2); end
if adddisp==1; for idsp=1:4; hold on; plot(kg,wgrdisp(:,idsp),'ko'); end; end

% Q
subplot(iy,ix,5);
term=xconQ; levels=levelsq;
if dolevs==1
contourf(kg,wg,term',levels, 'Linestyle', 'none'); 
caxis([levels(1) levels(end)]);% VERY IMPORTANT !!!!
colorbar('southoutside')
else
contourf(kg,wg,term','LineStyle','none');
end
xlim(xxrange); ylim(yyrange)
colorbar('southoutside')
xlabel('wavenumber(2$\pi$/40,000km)','interpreter','latex')
%  ylabel('frequency(cpd)','interpreter','latex')
title('(e) $q^{\prime}$','interpreter','latex');
if addlineintra==1; hold on; plot(wg*0,wg,'k-','linewidth',2); end
if addlineintra==1; hold on; plot(kg,kg*0+1/30,'k--','linewidth',2); end
if addlineintra==1; hold on; plot(kg,kg*0+1/90,'k--','linewidth',2); end
if adddisp==1; for idsp=1:4; hold on; plot(kg,wgrdisp(:,idsp),'ko'); end; end


end% dographpower
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART2: Power spectra TE, UZM

if dographspectrum==1
%  figure(31); clf; 
ix=3; iy=2;
subplot(iy,ix,6);
%
dologspectrumT=   1;% show log for the y axis
filetimesave=strcat('data_',datafolder,'/timeseries.nc');
xconT=ensospe_ncdfgetvar(filetimesave,'Te');% compute in ...timeserie_efficient
tcon=ensospe_ncdfgetvar(filetimesave,'ts')/365; % dimensionalize from day to years !!!!% compute in ...timeserie_efficient
%
varplot=xconT; % variable to plot
T=tcon(end)-tcon(1); % time of interval (days)
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
title('(f) $T_{E}$','interpreter','latex');
colorbar('location','southoutside','visible','off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Same power spectra for u average western half Pacific

subplot(iy,ix,3);

filetimesave=strcat('data_',datafolder,'/timeseries.nc');
xconT=ensospe_ncdfgetvar(filetimesave,'uw');% compute in ...timeserie_efficient
tcon=ensospe_ncdfgetvar(filetimesave,'ts')/365;% compute in ...timeserie_efficient
%
dologspectrumT=1; % do log spectrum in y

varplot=xconT-mean(xconT(:)); % variable to plot
T=tcon(end)-tcon(1); % time of interval (days)
N = length(tcon);% number of points
pow=abs(fft(varplot)); % absolute value of fft
pow=pow(1:N/2).^2; % take power of positive freq half
pow=pow./N^2; % normalise (note we should normalise by (N/2)^2 if cutting Nyquist)
freq=[0:N/2-1]/T; % find corresponding frequency in cpd (missing N/2 ?)
%  pow=log10(pow);



% 1) Cut to remove spurious values at very low w
%  wmincut=ensospe_searchclosest(freq,0.5*10e-5);
%  pow=pow(wmincut:end); freq=freq(wmincut:end);

if dologspectrumT==1
loglog(pow,freq);
else
plot(pow,freq);
end

% Compute extremely smoothed version
ntsmoo=     101; % t size of smooth (odd) (ref=25)
nsmoo=5;% number of times smoothing
for ismoo=1:nsmoo
pow=exp(ensospe_smoothn(log(pow),[1,ntsmoo]));% smooth
end
hold on; plot(pow,freq,'r-','linewidth',2);
%
yticks([0.01 0.1 1 10]);
xrange=[10e-11 1];
hold on; plot(xrange,xrange*0+1/30*365,'k--','linewidth',2);
hold on; plot(xrange,xrange*0+1/90*365,'k--','linewidth',2);
%  hold on; plot(xrange,xrange*0+1/10*365,'k-','linewidth',2);
%  hold on; plot(xrange,xrange*0+1/10,'k-','linewidth',2);

xlim(xrange)
if dologspectrumT==1
ylim([0.08,365/8]); %goes from 10 years to 5 days
else
ylim([0,15]); % intra
end
ylabel('$yr^{-1}$','interpreter','latex');
xlabel('$m^{2}.s^{-2}$','interpreter','latex');
%  ylabel('$Frequency (cpd)$','interpreter','latex');
title('(c) $u^{\prime}_{W}$','interpreter','latex');
colorbar('location','southoutside','visible','off');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end% dographspcetrum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
