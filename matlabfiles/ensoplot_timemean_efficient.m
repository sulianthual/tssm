% by Sulian Thual
% 
% plot timeseries from ENSO-MJO model solutions
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot parameters

%
ixmin=14; ixmax=28; % eastern Pacific half
ixshow=1;

dorecomputetimemean=        0           ;% recompute timemean
doshowtimemean =           1          ;% show timemean
doshowfluctuationstomean=        1;% add fluctuations mean on top


dolevs=        1           ; % do levels
levelsHa=    (0:0.1:1)*0.45-0.2  ; % (K.day-1)

filetimemeansave=strcat('data_',datafolder,'/timemean.nc');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recompute timemean

if dorecomputetimemean==1
% Read all files


% Index of each field
ixKa=1; ixRa=nx+1; ixA=nx*2+1; ixZ=nx*3+1;% index of variable start
ixKo=nx*4+1; ixRo=nx*4+nxo+1; ixT=nx*4+nxo*2+1; 
nxKa=nx; nxRa=nx*2; nxA=nx*3; nxZ=nx*4;% index of variable end
nxKo=nx*4+nxo; nxRo=nx*4+nxo*2; nxT=nx*4+nxo*3; 

% Loop  over restart files

Tmean=zeros(nxo,1); 
Kamean=zeros(nx,1);
Ramean=zeros(nx,1);
Zmean=zeros(nx,1); 
Amean=zeros(nx,1); 
MMAmean=zeros(nx,1); 
Komean=zeros(nxo,1);
Romean=zeros(nxo,1);

totalcount=(lastrestart-firstrestart+1); % total number of points
for indexrestart=firstrestart:lastrestart
strcat(['timemean,restart= ',num2str(indexrestart),' / ',num2str(lastrestart)])

% Read
dfolder=strcat('data_',datafolder); % setup is defined in main
fileout=strcat(dfolder,'/ensomjo_',num2str(indexrestart), '.nc');

XXs=ensospe_ncdfgetvar(fileout,'Xs');
XXs=phi0eq/totalcount*XXs;% divide by totalcount, take values at equator !
% Read and sum
term=XXs(ixT:nxT,:);
Tmean=Tmean+mean(term,2);
term=XXs(ixKa:nxKa,:);
Kamean=Kamean+mean(term,2);
term=XXs(ixRa:nxRa,:);
Ramean=Ramean+mean(term,2);
term=XXs(ixA:nxA,:);
Amean=Amean+mean(term,2);
term=XXs(ixZ:nxZ,:);
Zmean=Zmean+mean(term,2);
term=XXs(ixKo:nxKo,:);
Komean=Komean+mean(term,2);
term=XXs(ixRo:nxRo,:);
Romean=Romean+mean(term,2);

MMA=phi0eq/totalcount*ensospe_ncdfgetvar(fileout,'MMAs');
MMAmean=MMAmean+mean(MMA,2);

end

Hmean=Komean+Romean;
Umean=Komean-Romean;

term=so/phi0eq-HH*MMAmean;
kkkmean=ensonum_oldatm(3,nx,dx,10e-8,1, (1/2)*term); 
rrrmean=ensonum_oldatm(3,nx,dx,10e-8,1, -1*term);
MMumean=kkkmean-rrrmean; 
MMomean=-kkkmean-rrrmean; 
%  umean=ensonum_oldatm(3,nx,dx,10e-8,1, (3/2)*term); 
%  MMomean=ensonum_oldatm(3,nx,dx,10e-8,1, (3/2)*term); 

umean=(Kamean-Ramean);
omean=-(Kamean+Ramean);
qmean=Zmean-QQ*omean;

% write down
ensospe_ncdfmakevar(filetimemeansave,'T',{'xo','one'},Tmean,NaN,2); 
ensospe_ncdfmakevar(filetimemeansave,'H',{'xo','one'},Hmean,NaN,1); 
ensospe_ncdfmakevar(filetimemeansave,'U',{'xo','one'},Umean,NaN,1); 
ensospe_ncdfmakevar(filetimemeansave,'u',{'x','one'},umean,NaN,1); 
ensospe_ncdfmakevar(filetimemeansave,'a',{'x','one'},Amean,NaN,1); 
ensospe_ncdfmakevar(filetimemeansave,'theta',{'x','one'},omean,NaN,1);
ensospe_ncdfmakevar(filetimemeansave,'q',{'x','one'},qmean,NaN,1);
ensospe_ncdfmakevar(filetimemeansave,'abar',{'x','one'},MMAmean,NaN,1);
ensospe_ncdfmakevar(filetimemeansave,'ubar',{'x','one'},MMumean,NaN,1);
ensospe_ncdfmakevar(filetimemeansave,'thetabar',{'x','one'},MMomean,NaN,1);

end% recompute



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show timemean

if doshowtimemean==1

Tmean=ensospe_ncdfgetvar(filetimemeansave,'T');
Hmean=ensospe_ncdfgetvar(filetimemeansave,'H');
Umean=ensospe_ncdfgetvar(filetimemeansave,'U');
umean=ensospe_ncdfgetvar(filetimemeansave,'u');
omean=ensospe_ncdfgetvar(filetimemeansave,'theta');
qmean=ensospe_ncdfgetvar(filetimemeansave,'q');
Amean=ensospe_ncdfgetvar(filetimemeansave,'a');
MMAmean=ensospe_ncdfgetvar(filetimemeansave,'abar');
MMumean=ensospe_ncdfgetvar(filetimemeansave,'ubar');
MMomean=ensospe_ncdfgetvar(filetimemeansave,'thetabar');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Walker Circulation x-z
if 1==1

figure(iwind); iwind=iwind+1;, clf;

ix=2; iy=3; icount=1;
iwidd=1;

% ocean currents and winds (at equator)
%
% ubar
subplot(iy,ix,2); 
plot(xgo,ua*MMumean(1:28),'k-','linewidth',iwidd)
title('(b) Mean $\overline{u}$','interpreter','latex')
xlim([min(xgo) max(xgo)])
%  xlabel('x (1000km)','interpreter','latex');
ylabel('$m.s^{-1}$','interpreter','latex')
if doshowfluctuationstomean==1
hold on; plot(xgo,ua*MMumean(1:28)+ua*umean(1:28),'r--','linewidth',iwidd)
title('(b) Mean $\overline{u}+u^\prime$','interpreter','latex')
end


% U
if 0==1
subplot(iy,ix,5); 
plot(xgo,Ua*Umean,'k-','linewidth',iwidd)
title('(c) Mean $U $','interpreter','latex')
xlim([min(xgo) max(xgo)])
ylim([-0.5 0.5])
xlabel('x (1000km)','interpreter','latex');
ylabel('$m.s^{-1}$','interpreter','latex')
end

% H
subplot(iy,ix,4); 
plot(xgo,Ha*Hmean,'k-','linewidth',iwidd)
title('(c) Mean $H $','interpreter','latex')
xlim([min(xgo) max(xgo)])
ylim([-50 150])
%  xlabel('x (1000km)','interpreter','latex');
ylabel('$m$','interpreter','latex')

% Habar
subplot(iy,ix,4); 
plot(xgo,HH*oa/(ta/oneday)*MMAmean(1:28),'k-','linewidth',iwidd)
title('(b) Mean $\overline{H}\overline{a}$','interpreter','latex')
xlim([min(xgo) max(xgo)])
%  xlabel('x (1000km)','interpreter','latex');
ylabel('$m.s^{-1}$','interpreter','latex')
if doshowfluctuationstomean==1
hold on; plot(xgo,HH*oa/(ta/oneday)*MMAmean(1:28)+HH*oa/(ta/oneday)*Amean(1:28),'r--','linewidth',iwidd)
title('(b) Mean $\overline{H}(\overline{a}+a^{\prime})$','interpreter','latex')
end


% T
subplot(iy,ix,6); 
plot(xgo,Ta*Tmean,'k-','linewidth',iwidd)
title('(d) Mean $T$','interpreter','latex')
xlim([min(xgo) max(xgo)])

xlabel('x (1000km)','interpreter','latex');
ylabel('$K$','interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%55

% Compute Walker circulation sketch
zrange=[-1 16];
xrange=[min(xgo) max(xgo)];
zz=(0:0.05:1)*pi;
zzg=zz*16/pi;
MMuxmean=MMumean*0;
MMuxmean(1:nx-1)=MMumean(2:nx)-MMumean(1:nx-1); MMuxmean(nx)=MMumean(1)-MMumean(nx);
pass=xg(2)-xg(1); MMuxmean=-1/pass*MMuxmean;
MMuxmean=ua*16/pi/15000*MMuxmean;
afact=10000;% arbitrary fcator for arrows
MMuxmean=MMuxmean*afact; % amplification factor arbitrary


%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot walker circulation
subplot(iy,ix,[1,3,5]);
if dolevs==1; 
contourf(xg,zzg,HH*oa/(ta/oneday)*(MMAmean*sin(zz))',levelsHa,'linestyle','none')
caxis([levelsHa(1) levelsHa(end)]); 
%  caxis([0 levelsHa(end)]); 
else
contourf(xg,zzg,HH*oa/(ta/oneday)*(MMAmean*sin(zz))','linestyle','none')
end
%  contour(xg,zzg,HH*oa/(ta/oneday)*(MMAmean*sin(zz))','k-')
%  hold on; quiver(xg,zzg,(ua*MMumean*cos(zz))',(MMuxmean*sin(zz))','k-','autoscale','off','autoscalefactor',quivfact);
hold on; quiver(xg,zzg,(ua*MMumean*cos(zz))',(MMuxmean*sin(zz))','k-');
%  hold on; plot(xg,HH*oa/(ta/oneday)*MMAmean+HH*oa/(ta/oneday)*Amean,'r--')
%  hold on; plot([16.8750 16.8750],[-5 5],'k--')
xlabel('x (1000km)','interpreter','latex');
ylabel('z (km)','interpreter','latex');
xlim(xrange)
ylim(zrange)
title('(a) Mean ($\overline{H}\overline{a}$)','interpreter','latex');
hcb=colorbar('location','southoutside');
xlabel(hcb,'$K.day^{-1}$','interpreter','latex')

%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean state contours x-y
if 0==1

figure(iwind); iwind=iwind+1;, clf;

ix=3; iy=2; icount=1;

yrange=[-5 5];
yrangeo=[-1.5 1.5];
xrange=[min(xgo) max(xgo)];
% winds u
subplot(iy,ix,icount); icount=icount+1;
contourf(xg,yyg,(ua*MMumean*phi0/phi0eq)')
%  hold on; plot(xg,ua*MMumean+ua*umean),'r--')
%  hold on; plot([16.8750 16.8750],[-30 30],'k--')
xlabel('x (1000km)');
xlim(xrange)
ylim(yrange)
ylabel('y (1000km)');
title('$\overline{u}(m.s^{-1})$','interpreter','latex');
hcb=colorbar('location','southoutside');
%  xlabel(hcb,'$m.s^{-1}$', 'interpreter','latex'); 

% theta
subplot(iy,ix,icount); icount=icount+1;
contourf(xg,yyg,(oa*MMomean*phi0/phi0eq)')
%  hold on; plot(xg,oa*MMomean+oa*omean,'r--')
%  hold on; plot([16.8750 16.8750],[-30 30],'k--')
xlabel('x (1000km)');
xlim(xrange)
ylim(yrange)
title('$\overline{\theta}(K)$','interpreter','latex');
hcb=colorbar('location','southoutside');


% convection a
subplot(iy,ix,icount); icount=icount+1;
contourf(xg,yyg,HH*oa/(ta/oneday)*(MMAmean*phi0/phi0eq)')
%  hold on; plot(xg,HH*oa/(ta/oneday)*MMAmean+HH*oa/(ta/oneday)*Amean,'r--')
%  hold on; plot([16.8750 16.8750],[-5 5],'k--')
xlabel('x (1000km)');
xlim(xrange)
ylim(yrange)
title('$\overline{H}\overline{a}(K.day^{-1})$','interpreter','latex');
hcb=colorbar('location','southoutside');

% ocean currents
subplot(iy,ix,icount); icount=icount+1;
contourf(xgo,YYg,Ua*(Umean*psi0/phi0eq)')
xlabel('x (1000km)');
xlim(xrange)
ylim(yrangeo)
ylabel('y (1000km)');
title('$U(m.s^{-1})$','interpreter','latex');
hcb=colorbar('location','southoutside');

% thermocline
subplot(iy,ix,icount); icount=icount+1;
contourf(xgo,YYg,Ha*(Hmean*psi0/phi0eq)')
xlabel('x (1000km)');
xlim(xrange)
ylim(yrangeo)
title('$H(m)$','interpreter','latex');
hcb=colorbar('location','southoutside');

% sst
subplot(iy,ix,icount); icount=icount+1;
contourf(xgo,YYg,Ta*(Tmean*psi0/phi0eq)')

xlabel('x (1000km)');
xlim(xrange)
ylim(yrangeo)
title('$T(K)$','interpreter','latex');
hcb=colorbar('location','southoutside');


end


%
end% doshowtimemean
