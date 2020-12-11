% by Sulian Thual
% 
% contour hovmollers from ENSO-MJO model solutions
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME RANGE



indexrestart=1; nrestarts=2; % restart files span to use
trange=[NaN NaN];% timerange to use (Nan for all restart span)

% MAIN WP
if 1==1
indexrestart=401; nrestarts=2; trange=[1603 1608];%[1603 1608];%[485.75 487.25]*3.3;%doshowatmhov=1
indexrestart=404-10; nrestarts=21; trange=[480 500]*3.3;%[485.75 487.25];doshowatmhov=0
levels1proj=(-1.2:0.1:1.2)*0.1;%MJO
levels1projstd=(0:0.1:1)*0.06  ; % std(MJO)
end

% CRUDE INTRA
if 0==1
indexrestart=230; nrestarts=2; trange=[278.5 280]*3.3;%doshowatmhov=1
%  indexrestart=231-10; nrestarts=21; trange=[270 290]*3.3;%[NaN NaN];%doshowatmhov=0
levels1proj=(-1.2:0.1:1.2)*0.2; % MJO
levels1projstd=(0:0.1:1)*0.12  ; % std (MJO)
end

% MAIN WP new kde=8
if 1==1
indexrestart=51; nrestarts=2; trange=[207 212];%[NaN NaN];%[424.5 426]*3.3;%;%doshowatmhov=1
%  levels1proj=(-1.2:0.1:1.2)*0.1;%MJO
%  levels1projstd=(0:0.1:1)*0.06  ; % std(MJO)
end

% Just quick detection
if 1==1
%  indexrestart=178; nrestarts=2;  trange=[712.5 715]; % NOT TOO BAD wind bursts
%  indexrestart=234; nrestarts=2; trange=[937.5 939];% meh
%  indexrestart=251; nrestarts=1; trange=[1000.5 1002]; % not too bad
indexrestart=274; nrestarts=2; trange=[1095.6 1096.75];%[1095.5 1097]+[0.25 -0.25];% VERY GOOD WIND BURTS FOR THE DRAFT
%  trange=[NaN NaN];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% PARAMETERS 
hovremovetimemean=        1          ;% remove time mean in hovmoler (call timemean function before): Note Habar+Hamean adds the timemean again to ensure >0


showdrafthov=              1           ;% show either hovmoller for the draft
doshowatmhov=              1          ;% =1 draft hovmoller atm/intra oriented, =0 draft hovmoller inter atm
doshowatmhovtype3=          1; % doshowatmhov=1 and this=1 does 3rd type, showing ubar+u, habar+ha, ubar, habar, u,ha....
showHabarHahov =          1       ;% hos Habar+Ha instead (only for atm plot: FOR MAIN WP)
doshiftatmcircu=                1           ;% shift x-circularly atm to show Indian ocean
ndxcirc=20;% number of points shifted

showquickhovsst=           0       ;% just show sst quickly for event detection
showquickhovaddwinds=      1       ;% add uwinds for some nice wind bursts (longer)


xrangeatm=[0 16.8750]; % atm range same as ocean
xrangeatm=[0 40]; % atm range full equator
xrangeatm=[-10 16.8750]; % atm range including Indian ocean (need to shift circularly)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% regular contours

colormap(jet);

doensospe_filterkw=     0          ; % filter within k,w, for the hovmullers (ref=0) 
kmin=1; kmax=3; wmin=1/90; wmax=1/30; % in cpd, for hovmullers filtering


dolevs=        1           ; % do levels

levelsu=     (-1:0.1:1)*13  ; % (m.s-1) 
levelsubar=     (-1:0.1:1)*6  ; % (m.s-1) 
levelso=     (-1:0.1:1)*4+0.5  ; %  (K)
levelso=     (-1:0.1:1)*3 ; %  (K) CRUDE INTRA
levelsq=     (-1:0.1:1)*7  ; %(K)
levelsHa=    (-1:0.1:1)*1.2  ; % (K.day-1) CRUDE INTRA

levelsHabar=  (0:0.1:1)*0.3  ; % (K.day-1)
levelsU=     (-1:0.1:1)*0.75; % (m.s-1) 
levelsH=     (-1:0.1:1)*130; % (m) 
levelsT=     (-1:0.1:1)*7.5 ; % (K) 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Index of each field
ixKa=1; ixRa=nx+1; ixA=nx*2+1; ixZ=nx*3+1;% index of variable start
ixKo=nx*4+1; ixRo=nx*4+nxo+1; ixT=nx*4+nxo*2+1; 
nxKa=nx; nxRa=nx*2; nxA=nx*3; nxZ=nx*4;% index of variable end
nxKo=nx*4+nxo; nxRo=nx*4+nxo*2; nxT=nx*4+nxo*3; 

% Read

dfolder=strcat('data_',datafolder); % setup is defined in main
fileout=strcat(dfolder,'/ensomjo_',num2str(indexrestart), '.nc');
ts=ensospe_ncdfgetvar(fileout,'ts');
XXs=ensospe_ncdfgetvar(fileout,'Xs');
MMAs=ensospe_ncdfgetvar(fileout,'MMAs');

% concatenate
if nrestarts>1
for irestart=indexrestart+1:indexrestart+nrestarts
dfolder=strcat('data_',datafolder); % setup is defined in main
fileout=strcat(dfolder,'/ensomjo_',num2str(irestart), '.nc');
ts2=ensospe_ncdfgetvar(fileout,'ts');
XXs2=ensospe_ncdfgetvar(fileout,'Xs');
MMA2s=ensospe_ncdfgetvar(fileout,'MMAs');

XXs=[XXs';XXs2']';
MMAs=[MMAs';MMA2s']';
ts=[ts, ts2];
end
end

nts=length(ts);% change total time
ts=ts*3.3; % Dimensionalizes time to days !!!!! (timestep was each 3.3 days)




% Attribute
Kas=XXs(ixKa:nxKa,:);
Ras=XXs(ixRa:nxRa,:);
As=XXs(ixA:nxA,:);
Zs=XXs(ixZ:nxZ,:);
Kos=XXs(ixKo:nxKo,:);
Ros=XXs(ixRo:nxRo,:);
Ts=XXs(ixT:nxT,:);

% reconstruct winds
mus=zeros(nx,nts);
for its=1:nts
Eq=zeros(nx,1); Eq(1:nxo,1)=alphaq*Ts(:,its); 
Eqzma=Eq-mean(Eq);
Abar=1/(HH*(1-QQ))*(chia*Eqzma+sq/phi0eq-QQ*so/phi0eq);
Abar(Abar<=0)=10e-5;% impose Abar>0
term=so/phi0eq-HH*Abar;
mus(:,its)=ensonum_oldatm(3,nx,dx,10e-8,1, (3/2)*term); 
end

% Reconstruct variables at equator
us=(Kas-Ras)*phi0eq;
os=-(Kas+Ras)*phi0eq;
aas=As*phi0eq;
qs=QQ*os+Zs*phi0eq;
Ts=Ts*phi0eq;
Hs=(Kos+Ros)*psi0eq;
Us=(Kos-Ros)*psi0eq;
mus=mus*psi0eq;
MMAs=MMAs*phi0eq;

if hovremovetimemean==1
filetimemeansave=strcat('data_',datafolder,'/timemean.nc');
Tmean=ensospe_ncdfgetvar(filetimemeansave,'T');
Hmean=ensospe_ncdfgetvar(filetimemeansave,'H');
Amean=ensospe_ncdfgetvar(filetimemeansave,'a');
Umean=ensospe_ncdfgetvar(filetimemeansave,'U');
umean=ensospe_ncdfgetvar(filetimemeansave,'u');
omean=ensospe_ncdfgetvar(filetimemeansave,'theta');
qmean=ensospe_ncdfgetvar(filetimemeansave,'q');
MMAmean=ensospe_ncdfgetvar(filetimemeansave,'abar');
MMumean=ensospe_ncdfgetvar(filetimemeansave,'ubar');
%  MMomean=ensospe_ncdfgetvar(filetimemeansave,'thetabar');
us=us-repmat(umean,[1,nts]);
os=os-repmat(omean,[1,nts]);
qs=qs-repmat(qmean,[1,nts]);
aas=aas-repmat(Amean,[1,nts]);
Ts=Ts-repmat(Tmean,[1,nts]);
Hs=Hs-repmat(Hmean,[1,nts]);
Us=Us-repmat(Umean,[1,nts]);
MMAs=MMAs-repmat(MMAmean,[1,nts]);
mus=mus-repmat(MMumean,[1,nts]);

end




% Filter kw in atmosphere
if doensospe_filterkw==1

for ivarfilt=1:4
if ivarfilt==1; term=us; end
if ivarfilt==2; term=os; end
if ivarfilt==3; term=qs; end
if ivarfilt==4; term=aas; end
pass=size(term); nth=pass(2);
term=fftshift(fft2(term)); 
kg=ensospe_fftkspe(nx,dx); % adim
wg=ensospe_fftkspe(nth,dt*mts); %adim
kg=kg/xa*40000*1000/(2*pi);
wg=wg/ta*oneday; 
wgmin=ensospe_searchclosest(wg,wmin); wgmax=ensospe_searchclosest(wg,wmax);
kgmin=ensospe_searchclosest(kg,-kmax); kgmax=ensospe_searchclosest(kg,-kmin);% BEWARE AXIS IS BACKWARD
term2=term*0; term2(kgmin:kgmax,wgmin:wgmax)=term(kgmin:kgmax,wgmin:wgmax); clear term;
term2=ifft2(ifftshift(term2)); term2=real(term2);
if ivarfilt==1; us=term2; end 
if ivarfilt==2; os=term2; end 
if ivarfilt==3; qs=term2; end 
if ivarfilt==4; aas=term2; end 
clear term2;
end

end



%
ts=ts/365;% years

% Apply instead a running std
if 0==1
windstd=45;
[pass,ntoner]=size(us);
usstd=us*0;
for ix=1:nx
for it=1+windstd:ntoner-windstd
usstd(ix,it)=std(us(ix,it-windstd:it+windstd));
end
end
us=usstd; 
clear usstd
end

% range time (default if Nan chosen)
if isfinite(trange(1))==0
trange=[min(ts) max(ts)];
end
trangeocean=trange;

% shift circularly atmosphere variables
if doshiftatmcircu==1
us=circshift(us,ndxcirc,1);
os=circshift(os,ndxcirc,1);
qs=circshift(qs,ndxcirc,1);
aas=circshift(aas,ndxcirc,1);
Kas=circshift(Kas,ndxcirc,1);
Ras=circshift(Ras,ndxcirc,1);
Zs=circshift(Zs,ndxcirc,1);
As=circshift(As,ndxcirc,1);
mus=circshift(mus,ndxcirc,1);
MMAs=circshift(MMAs,ndxcirc,1);
xgsh=circshift(xg,ndxcirc);
xgsh(1:ndxcirc)=xgsh(1:ndxcirc)-40;
if hovremovetimemean==1% should circle all means, lazy just do the needed ones
Amean=circshift(Amean,ndxcirc,1);
MMAmean=circshift(MMAmean,ndxcirc,1);
end
%  stop
else
xgsh=xg;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Hovmollers of variables: type=1

if showdrafthov==1

figure(iwind); iwind=iwind+1;, clf;
colormap(jet);


%%%%%%%%%%%%%%%% u
if doshowatmhov; ix=8; iy=1; icount=1; else ix=6; iy=1; icount=1; end
if doshowatmhov==1
% u
subplot(iy,ix,2); icount=icount+1; 
term=ua*us; levels=levelsu; titlepass='(b) $u^{\prime}$'; 
if dolevs; contourf(xgsh,ts,term',levels,'LineStyle','none');
caxis([levels(1) levels(end)]);% VERY IMPORTANT !!!!
else contourf(xgsh,ts,term','LineStyle','none'); end
ylim(trange); 
xlim(xrangeatm);
xlabel('x(1000km)','interpreter','latex')
hold on; plot([0 0],trange,'r-','linewidth',1);
%  ylabel('time (days)')
%  ylabel('time (years)')
%  yticks([]);
yticklabels([' '])
title(titlepass,'interpreter','latex')
hcb=colorbar('location','southoutside');
xlabel(hcb,'$m.s^{-1}$', 'interpreter','latex'); 
end
%
%%%%%%%%%%%%%%%% theta
if doshowatmhov==1
% theta
subplot(iy,ix,3); icount=icount+1; 
term=oa*os; levels=levelso; titlepass='(c) $\theta^{\prime}$'; 
if dolevs; contourf(xgsh,ts,term',levels,'LineStyle','none'); 
caxis([levels(1) levels(end)]);
else contourf(xgsh,ts,term','LineStyle','none'); end
if doshowatmhovtype3==1; hold on; contour(xgsh,ts,term',[0,0],'k--'); end
ylim(trange); 
xlim(xrangeatm);
xlabel('x(1000km)','interpreter','latex')
hold on; plot([0 0],trange,'r-','linewidth',1);
yticklabels([' '])
%  ylabel('time (days)')
title(titlepass,'interpreter','latex')
hcb=colorbar('location','southoutside');
xlabel(hcb,'$K$', 'interpreter','latex'); 
end

%%%%%%%%%%%%%%%% type3: plot ubar+u (instead of theta)
if doshowatmhov==1
if doshowatmhovtype3==1
% theta
subplot(iy,ix,3);
term=ua*us+ua*mus; levels=levelsu; titlepass='(c) $\overline{u}+u^{\prime}$'; 
if dolevs; contourf(xgsh,ts,term',levels,'LineStyle','none');
caxis([levels(1) levels(end)]);% VERY IMPORTANT !!!!
else contourf(xgsh,ts,term','LineStyle','none'); end
if doshowatmhovtype3==1; hold on; contour(xgsh,ts,term',[0,0],'k--'); end
ylim(trange); 
xlim(xrangeatm);
xlabel('x(1000km)','interpreter','latex')
hold on; plot([0 0],trange,'r-','linewidth',1);
%  ylabel('time (days)')
%  ylabel('time (years)')
%  yticks([]);
yticklabels([' '])
title(titlepass,'interpreter','latex')
hcb=colorbar('location','southoutside');
xlabel(hcb,'$m.s^{-1}$', 'interpreter','latex'); 
end
end

%%%%%%%%%%%%%%%% q
if doshowatmhov==1
% qs
subplot(iy,ix,4); icount=icount+1; 
term=oa*qs; levels=levelsq; titlepass='(d) $q^{\prime}$'; 
if dolevs; contourf(xgsh,ts,term',levels,'LineStyle','none'); 
caxis([levels(1) levels(end)]);
else contourf(xgsh,ts,term','LineStyle','none'); end
ylim(trange); 
xlim(xrangeatm);
xlabel('x(1000km)','interpreter','latex')
hold on; plot([0 0],trange,'r-','linewidth',1);
yticklabels([' '])
%  ylabel('time (days)')
title(titlepass,'interpreter','latex')
hcb=colorbar('location','southoutside');
xlabel(hcb,'$K$', 'interpreter','latex'); 
end
%

%%%%%%%%%%%%%%%% Ha (or HA+Habar) instead of q
if doshowatmhov==1
if doshowatmhovtype3==1
subplot(iy,ix,4); 

term=HH*oa/(ta/oneday)*(aas+MMAs); 
if hovremovetimemean==1;% reput the timemean here if taken out
term=term+HH*oa/(ta/oneday)*(repmat(Amean,[1,nts])+repmat(MMAmean,[1,nts]));
end
% some tuning on colorbar
levels=(0:0.1:1)    ;%[(-0.5:0.1:1),2]  ; % (K.day-1) TOTAL Ha+Habar, WP MAIN
levels2=(0:0.1:1);
amap=colormap(jet);
amap=amap(16:end,:); % restrict colormap
colormap(gca,amap);% change only this subplot
titlepass='(e) $\overline{H}(a^{\prime}+\overline{a})$'; 

if dolevs; contourf(xgsh,ts,term',levels,'LineStyle','none'); 
caxis([levels2(1) levels2(end)]);
else contourf(xgsh,ts,term','LineStyle','none'); end
ylim(trange); 
xlim(xrangeatm);
xlabel('x(1000km)','interpreter','latex')
hold on; plot([0 0],trange,'r-','linewidth',1);
yticklabels([' '])
%  ylabel('time (days)')
title(titlepass,'interpreter','latex')
hcb=colorbar('location','southoutside');
xlabel(hcb,'$K.day^{-1}$', 'interpreter','latex'); 
%  colormap(jet);% reput normal
end
end


%%%%%%%%%%%%%%%% Ha (or HA+Habar)
if doshowatmhov==1
subplot(iy,ix,5); icount=icount+1; 

% Show Habar+Ha or Habar only
if showHabarHahov==1
%
term=HH*oa/(ta/oneday)*(aas+MMAs); 
if hovremovetimemean==1;% reput the timemean here if taken out
term=term+HH*oa/(ta/oneday)*(repmat(Amean,[1,nts])+repmat(MMAmean,[1,nts]));
end
% remove value below 0.1
%  term(term<0.1)=NaN;
% some tuning on colorbar
levels=(0:0.1:1)    ;%[(-0.5:0.1:1),2]  ; % (K.day-1) TOTAL Ha+Habar, WP MAIN
levels2=(0:0.1:1);
amap=colormap(jet);
amap=amap(16:end,:); % restrict colormap
colormap(gca,amap);% change only this subplot
titlepass='(e) $\overline{H}(a^{\prime}+\overline{a})$'; 
else
% Has
term=HH*oa/(ta/oneday)*aas; levels=levelsHa; titlepass='(e) $\overline{H}a^{\prime}$'; 
levels2=levels;
end
if dolevs; contourf(xgsh,ts,term',levels,'LineStyle','none'); 
caxis([levels2(1) levels2(end)]);
else contourf(xgsh,ts,term','LineStyle','none'); end
ylim(trange); 
xlim(xrangeatm);
xlabel('x(1000km)','interpreter','latex')
hold on; plot([0 0],trange,'r-','linewidth',1);
yticklabels([' '])
%  ylabel('time (days)')
title(titlepass,'interpreter','latex')
hcb=colorbar('location','southoutside');
xlabel(hcb,'$K.day^{-1}$', 'interpreter','latex'); 
%  colormap(jet);% reput normal
end



%%%%%%%%%%%%%%%% ubar
if doshowatmhov==1
subplot(iy,ix,6); icount=icount+1; titlepass='(f) $\overline{u}$'; 
else
% ubar
subplot(iy,ix,3); icount=icount+1; titlepass='(c) $\overline{u}$'; 
end

term=ua*mus; levels=levelsu; 
if dolevs; contourf(xgsh,ts,term',levelsubar,'LineStyle','none');
caxis([levelsubar(1) levelsubar(end)]);% VERY IMPORTANT !!!!
else contourf(xgsh,ts,term','LineStyle','none'); end
if doshowatmhovtype3==1; hold on; contour(xgsh,ts,term',[0,0],'k--'); end
ylim(trange); 
xlim([0 16.8750]);
xlabel('x(1000km)','interpreter','latex')
yticklabels([' '])
%  ylabel('time (days)')
%  ylabel('time (years)')
title(titlepass,'interpreter','latex')
hcb=colorbar('location','southoutside');
xlabel(hcb,'$m.s^{-1}$', 'interpreter','latex'); 


%%%%%%%%%%%%%%%% Habar
if doshowatmhov==1
else
% Habars
subplot(iy,ix,2); icount=icount+1; 
%  MMAs(MMAs<=0)=10e-5;
term=HH*oa/(ta/oneday)*MMAs; levels=levelsHabar; titlepass='(b) $\overline{H}\overline{a}$'; 
%  term=log(term);% TEST LOG SCALE
if dolevs; contourf(xgsh,ts,term',levels,'LineStyle','none'); 
caxis([levels(1) levels(end)]);
else contourf(xgsh,ts,term','LineStyle','none'); end
ylim(trange); 
xlim([0 16.8750]);
yticklabels([' '])
xlabel('x(1000km)','interpreter','latex')
%  ylabel('time (years)')
title(titlepass,'interpreter','latex')
hcb=colorbar('location','southoutside');
xlabel(hcb,'$K.day^{-1}$', 'interpreter','latex'); 
end

%%%%%%%%%%%%%%%% U
if 0==1
% U
subplot(iy,ix,4); icount=icount+1; 
term=Ua*Us; levels=levelsU; titlepass='(d) $U$'; 
if dolevs; contourf(xgo,ts,term',levels,'LineStyle','none'); 
caxis([levels(1) levels(end)]);
else contourf(xgo,ts,term','LineStyle','none'); end
ylim(trangeocean); 

yticklabels([' '])
%  ylabel('time (years)')
xlabel('x(1000km)','interpreter','latex')
%  ylabel('time (days)')
title(titlepass,'interpreter','latex')
hcb=colorbar('location','southoutside');
xlabel(hcb,'$m.s^{-1}$', 'interpreter','latex'); 
end

%%%%%%%%%%%%%%%% H
if doshowatmhov==1
subplot(iy,ix,7); icount=icount+1; titlepass='(g) $H$';
else
subplot(iy,ix,4); icount=icount+1; titlepass='(d) $H$';
end
% H
term=Ha*Hs; levels=levelsH;  
if dolevs; contourf(xgo,ts,term',levels,'LineStyle','none'); 
caxis([levels(1) levels(end)]);
else contourf(xgo,ts,term','LineStyle','none'); end
ylim(trangeocean);
yticklabels([' '])
xlabel('x(1000km)','interpreter','latex')
%  ylabel('time (days)')
title(titlepass,'interpreter','latex')
hcb=colorbar('location','southoutside');
xlabel(hcb,'$m$', 'interpreter','latex'); 

%%%%%%%%%%%%%%%% T
if doshowatmhov==1
subplot(iy,ix,8); icount=icount+1; titlepass='(h) $T$'; 
else
subplot(iy,ix,5); icount=icount+1; titlepass='(e) $T$'; 
end
% T
subplot(iy,ix,icount); icount=icount+1; 
term=Ta*Ts; levels=levelsT; 
if dolevs; contourf(xgo,ts,term',levels,'LineStyle','none'); 
caxis([levels(1) levels(end)]);
else contourf(xgo,ts,term','LineStyle','none'); end
ylim(trangeocean); 
yticklabels([' '])
xlabel('x(1000km)','interpreter','latex')
%  ylabel('time (days)')
title(titlepass,'interpreter','latex')
hcb=colorbar('location','southoutside');
xlabel(hcb,'$K$', 'interpreter','latex'); 


% Part II: timeserie of TE
if doshowatmhov==1

else
subplot(iy,ix,6); 
filetimesave=strcat('data_',datafolder,'/timeseries.nc');
xconT=ensospe_ncdfgetvar(filetimesave,'Te');% compute in ...timeserie_efficient
xconTw=ensospe_ncdfgetvar(filetimesave,'Tw');% compute in ...timeserie_efficient
tcon=ensospe_ncdfgetvar(filetimesave,'ts');% compute in ...timeserie_efficient
plot(xconT,tcon/365,'k-','linewidth',1); 
hold on; plot(xconTw,tcon/365,'r-','linewidth',1); 
yticklabels([' '])
ylim(trange);
xlim([-9 9]);
xlabel('K','interpreter','latex');
title('(f) $T_{E},T_{W}$','interpreter','latex');
colorbar('location','southoutside','visible','off');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data projection


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data projection on mode
%
imodshow= 2; % eignemode from linear stability on which to project 

%
doensospe_filterkw2=     1          ; % filter within k-w the signal before projection (ref=1)
kmin2=1; kmax2=3;% in (2pi/40000km)
wmin2=1/90; wmax2=1/30; %
%
donullmean=     0          ; % put k=0 to zero when projecting (ref=0) 
donullnyquist=  0          ; % put nyquist k to zero (ik=1) when projecting (ref=0)
dorefphase=     1          ; % impose reference phase (ref=1)
%
dolevsproj=       1       ; % levels for data projection
%  levels1proj=(-1.2:0.1:1.2)*0.1  ; % real, see in main for def
%  levels1proj=(-1.2:0.1:1.2)*0.2  ; % real, GOOD CRUDE INTRA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read Outputs 
%  Kas=XXs(ixKa:nxKa,:);
%  Ras=XXs(ixRa:nxRa,:);
%  As=XXs(ixA:nxA,:);
%  Zs=XXs(ixZ:nxZ,:);

pass=size(Kas); nth=pass(2);

% Filter kw in atmosphere
if doensospe_filterkw2==1

for ivarfilt=1:4
if ivarfilt==1; term=Kas; end
if ivarfilt==2; term=Ras; end
if ivarfilt==3; term=Zs; end
if ivarfilt==4; term=As; end
pass=size(term); nth=pass(2);
term=fftshift(fft2(term)); 
kg=ensospe_fftkspe(nx,dx); % adim
wg=ensospe_fftkspe(nth,dt*mts); %adim
kg=kg/xa*40000*1000/(2*pi);
wg=wg/ta*oneday; 
wgmin=ensospe_searchclosest(wg,wmin); wgmax=ensospe_searchclosest(wg,wmax);
kgmin=ensospe_searchclosest(kg,-kmax); kgmax=ensospe_searchclosest(kg,-kmin);% BEWARE AXIS IS BACKWARD
term2=term*0; term2(kgmin:kgmax,wgmin:wgmax)=term(kgmin:kgmax,wgmin:wgmax); clear term;
term2=ifft2(ifftshift(term2)); term2=real(term2);
if ivarfilt==1; Kas=term2; end 
if ivarfilt==2; Ras=term2; end 
if ivarfilt==3; Zs=term2; end 
if ivarfilt==4; As=term2; end 
clear term2;
end

end

% Compute FFT zonal

for ivar=1:4
if ivar==1; x=Kas; end
if ivar==2; x=Ras; end
if ivar==3; x=Zs; end
if ivar==4; x=As; end
for kts=1:nth; x(:,kts)=fftshift(fft(x(:,kts))); end; % zonal
if donullmean==1; x(nx/2+1,:)=0; end
if donullnyquist==1; x(1,:)=0; end
if ivar==1; Kas=x; end
if ivar==1; Ras=x; end
if ivar==1; Zs=x; end
if ivar==1; As=x; end
end

%
% Compute stability
kk=ensospe_fftkspe(nx,dx); 
xproj=zeros(nx,nth);

for ik=1:length(kk)
ekr=0.1*(1/(30*oneday))*ta; % assumed dissip all variables, must have e=0.1 cf dim time=33 days
sorefdisp=2.2; % assumed sq/so here

% compute stability
[wr,wi,vects]=ensospe_stabskelnew(kk(ik),0.1*(1/(30*oneday))*ta,QQ,GG,HH,sorefdisp);
wr=wr(imodshow); 
wrtest(ik)=wr;
vects=squeeze(vects(:,imodshow)); % X=K,Rm,Am,Zm
%
% Modifiy phase (convention here is that phase(K)=0)
if dorefphase==1
term=vects(1); phaseK=angle(term);
vects=vects.*exp(-ii*phaseK);
end
%
% compute data projection
for kts=1:nth
% vector from simu (with zonal FFT)
Kp=Kas(ik,kts); Rmp=Ras(ik,kts); Ap=As(ik,kts); Zp=Zs(ik,kts); 
nxs=4;
Xsim=[Kp,Rmp,Ap,Zp].'; 
xproj(ik,kts)=Xsim.'*vects; 
end% kts
end% ik
%
% Inverse Fourier
for kts=1:nth; xproj(:,kts)=ifft(ifftshift(xproj(:,kts))); end; 


% Graph data projection 
if doshowatmhov==1

%figure(iwind); iwind=iwind+1;, clf;
subplot(iy,ix,1); icount=icount+1;
term=real(xproj);
if dolevsproj; contourf(xgsh,ts,term',levels1proj,'LineStyle','none');
caxis([levels1proj(1) levels1proj(end)]);
else contourf(xgsh,ts,term','LineStyle','none'); end

colorbar('location','southoutside');
xlabel('x(1000km)','interpreter','latex')
hold on; plot([0 0],trange,'r-','linewidth',1);
ylim(trange);
xlim(xrangeatm);
%  set(gca,'YTick',[]);
ylabel('time (years)','interpreter','latex');
title('(a) MJO','interpreter','latex')
end



% Graph data projection (absolute, smoothed)
if doshowatmhov==0
subplot(iy,ix,1); icount=icount+1;
term=real(xproj); 




% smoothing of absolute value
if 0==1
nsmoo=5;
nxsmoo=1; ntsmoo=101;
term=abs(xproj);
for ismoo=1:nsmoo; term=ensospe_smoothn(term,[nxsmoo,ntsmoo]); end
end


% Compute running standard deviation instead
if 1==1
term2=term*0;
nwind=450;
for it=1+nwind:length(ts)-nwind; 
for ix=1:nx
term2(ix,it)=std(term(ix,it-nwind:it+nwind));
end
end
term=term2; clear term2;
end

if dolevsproj; contourf(xgsh,ts,term',levels1projstd,'LineStyle','none');
caxis([levels1projstd(1) levels1projstd(end)]);
else contourf(xgsh,ts,term','LineStyle','none'); end

colorbar('location','southoutside');
xlabel('x(1000km)','interpreter','latex')
hold on; plot([0 0],trange,'r-','linewidth',1);
ylim(trange);
xlim(xrangeatm);
%  set(gca,'YTick',[]);
ylabel('time (years)','interpreter','latex');
title('(a) std(MJO)','interpreter','latex')


end

end% do draft graph


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% instead just show sst quickly to detect events


if showquickhovsst==1


figure(43);  clf;
if showquickhovaddwinds==1; subplot(1,2,1); else; end
colormap(jet);

titlepass='$T$'; 
% T
term=Ta*Ts; levels=levelsT; 
if dolevs; contourf(xgo,ts,term',levels,'LineStyle','none'); 
caxis([levels(1) levels(end)]);
else contourf(xgo,ts,term','LineStyle','none'); end
ylim(trangeocean); 
%  yticklabels([' '])
xlabel('x(1000km)','interpreter','latex')
%  ylabel('time (days)')
title(titlepass,'interpreter','latex')
hcb=colorbar('location','southoutside');
xlabel(hcb,'$K$', 'interpreter','latex'); 

if showquickhovaddwinds==1;% add winds
subplot(1,2,2)
term=ua*us; levels=levelsu; titlepass='(b) $u^{\prime}$'; 
if dolevs; contourf(xgsh,ts,term',levels,'LineStyle','none');
caxis([levels(1) levels(end)]);% VERY IMPORTANT !!!!
else contourf(xgsh,ts,term','LineStyle','none'); end
ylim(trange); 
xlim(xrangeatm);
xlabel('x(1000km)','interpreter','latex')
hold on; plot([0 0],trange,'r-','linewidth',1);

%  yticklabels([' '])
title(titlepass,'interpreter','latex')
hcb=colorbar('location','southoutside');
xlabel(hcb,'$m.s^{-1}$', 'interpreter','latex'); 
end

end
