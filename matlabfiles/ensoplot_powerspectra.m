% by Sulian Thual
% 
% contour power spectra in the atmosphere
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%
dorecomputepower==1;% RECOMPUTE POWER (VERY LONG !!!!)
%
if dorecomputepower==1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read all files

% Index of each field
ixKa=1; ixRa=nx+1; ixA=nx*2+1; ixZ=nx*3+1;% index of variable start
ixKo=nx*4+1; ixRo=nx*4+nxo+1; ixT=nx*4+nxo*2+1; 
nxKa=nx; nxRa=nx*2; nxA=nx*3; nxZ=nx*4;% index of variable end
nxKo=nx*4+nxo; nxRo=nx*4+nxo*2; nxT=nx*4+nxo*3; 

% Loop  over restart files

xconKa=((1:nx)*0)';
xconRa=((1:nx)*0)';
xconA=((1:nx)*0)';
xconZ=((1:nx)*0)';

for indexrestart=firstrestart+1:lastrestart
strcat(['restart= ',num2str(indexrestart),' / ',num2str(lastrestart)])

% Read
dfolder=strcat('data_',datafolder); % setup is defined in main
fileout=strcat(dfolder,'/ensomjo_',num2str(indexrestart), '.nc');


% Concatenate variables
XXs=ensospe_ncdfgetvar(fileout,'Xs');
Ka=XXs(ixKa:nxKa,:);
Ra=XXs(ixRa:nxRa,:);
A=XXs(ixA:nxA,:);
Z=XXs(ixZ:nxZ,:);
%  Ko=XXs(ixKo:nxKo,:);
%  Ro=XXs(ixRo:nxRo,:);
%  T=XXs(ixT:nxT,:);

xconKa=[xconKa,Ka];  
xconRa=[xconRa,Ra]; 
xconA=[xconA,A]; 
xconZ=[xconZ,A]; 

end

% Remove first dummy
xconKa=xconKa(:,2:end); 
xconRa=xconRa(:,2:end); 
xconA=xconA(:,2:end); 
xconZ=xconZ(:,2:end); 

xconu=xconKa-xconRa; 
xcono=-xconKa-xconRa; 
xconQ=xconZ+QQ*(xconKa+xconRa);
clear xconZ

% remove temporal mean
for ix=1:nx
xconu(ix,:)=xconu(ix,:)-mean(xconu(ix,:));
xcono(ix,:)=xcono(ix,:)-mean(xcono(ix,:));
xconA(ix,:)=xconA(ix,:)-mean(xconA(ix,:));
xconQ(ix,:)=xconQ(ix,:)-mean(xconQ(ix,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute power




pass=size(xconKa); nth=pass(2);
for ivar=1:4
if ivar==1; term=xconu; end
if ivar==2; term=xcono; end
if ivar==3; term=xconQ; end
if ivar==4; term=xconA; end

term=fftshift(fft2(term)); 
term=abs(term).^2/(nx*nth)^2; 
kg=ensospe_fftkspe(nx,dx); % adim
wg=ensospe_fftkspe(nth,dt*mts); %adim
kg=kg/xa*40000*1000/(2*pi);
wg=wg/ta*oneday; % BEWARE ----- DIRTY FIX
kg=-kg; % why ? this is no longer due to transpositions '(instead of .') of complex matrix ?
term=log10(term);
%
% Further treatment: smooth
nxsmoo=     1; % x size of smooth (odd) (ref=1)
ntsmoo=     25; % t size of smooth (odd) (ref=25)
term=ensospe_smoothn(term,[nxsmoo,ntsmoo]);% smooth
%
% Further treatment: undersample
nku=1; % k-space between samples 
nwu=11; % w-space between samples 
term=term(1:nku:end,1:nwu:end); 
kg=kg(1:nku:end);
wg=wg(1:nwu:end); 
%
if ivar==1; xconu=term; end
if ivar==2; xcono=term; end
if ivar==3; xconQ=term; end
if ivar==4; xconA=term; end
end

end% REDO COMPUTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contour

% restrict range
xxrange=[-5,5]; yyrange=[0,0.1];
%  xxrange=[-10,10]; yyrange=[0,0.2];

%  kgmin=ensospe_searchclosest(kg,-5); kgmax=ensospe_searchclosest(kg,5);
%  wgmin=ensospe_searchclosest(wg,0); wgmax=ensospe_searchclosest(wg,0.1);
%  kgmin=28;kgmax=38;
%  wgmin=2191; wgmax=3651;
%  kg=kg(kgmin:kgmax); wg=wg(wgmin:wgmax);
%  xconu=xconu(kgmin:kgmax,wgmin:wgmax);
%  xcono=xconu(kgmin:kgmax,wgmin:wgmax);
%  xconQ=xconu(kgmin:kgmax,wgmin:wgmax);
%  xconA=xconu(kgmin:kgmax,wgmin:wgmax);

%%%%%%%%%%%%%%%%%
figure(iwind); clf; iwind=iwind+1;
ix=2; iy=2; icount=1;
colormap(jet)
% u
subplot(iy,ix,icount); icount=icount+1;
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
xlabel('wavenumber(2pi/40,000km)')
ylabel('frequency(cpd)')
title('$u$','interpreter','latex');

% theta
subplot(iy,ix,icount); icount=icount+1;
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
xlabel('wavenumber(2pi/40,000km)')
ylabel('frequency(cpd)')
title('$\theta$','interpreter','latex');

% Q
subplot(iy,ix,icount); icount=icount+1;
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
xlabel('wavenumber(2pi/40,000km)')
ylabel('frequency(cpd)')
title('$q$','interpreter','latex');

% A
subplot(iy,ix,icount); icount=icount+1;
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
xlabel('wavenumber(2pi/40,000km)')
ylabel('frequency(cpd)')
title('$a$','interpreter','latex');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0==1

xxrange=[-5,5]; yyrange=[0,0.05];

figure(iwind); clf; iwind=iwind+1;
term=xconA; levels=-9:0.1:-6;
if 1==1
contourf(kg,wg,term',levels, 'Linestyle', 'none'); 
caxis([levels(1) levels(end)]);% VERY IMPORTANT !!!!
colorbar('southoutside')
else
contourf(kg,wg,term','LineStyle','none');
end
hold on; plot(kg,kg*0+0.01,'k--');
xlim(xxrange); ylim(yyrange)
colorbar('southoutside')
xlabel('wavenumber(2pi/40,000km)')
ylabel('frequency(cpd)')
title('$u_{SDE}$','interpreter','latex');


figure(iwind); clf; iwind=iwind+1;
term=xconu; levels=levelsa;
if 1==1
contourf(kg,wg,term',levels, 'Linestyle', 'none'); 
caxis([levels(1) levels(end)]);% VERY IMPORTANT !!!!
colorbar('southoutside')
else
contourf(kg,wg,term','LineStyle','none');
end
hold on; plot(kg,kg*0+0.01,'k--');
xlim(xxrange); ylim(yyrange)
colorbar('southoutside')
xlabel('wavenumber(2pi/40,000km)')
ylabel('frequency(cpd)')
title('$u$','interpreter','latex');


end

