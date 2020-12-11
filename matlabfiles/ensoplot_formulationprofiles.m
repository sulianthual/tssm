% Model Formulation: plot some profiles
% by Sulian Thual

figure(32); clf;
ix=2; iy=2; 

% Plot profiles sq, so: Warm Pool simulation
if 1==1
subplot(iy,ix,1);
sq=          2.2*(1+0.6*cos(2*pi/LL*xx)'); % warm pool
so=          sq; % so in skeleton model at equator
plot(xg,oa/ta/0.1*oneday*sq,'k-'); hold on;
plot(xg,oa/ta/0.1*oneday*so,'r--'); hold on;
plot([xg(28) xg(28)],[0,2],'k--','linewidth',2);
%  plot(xg,(sq-QQ*so)/(1-QQ),'b-'); hold on;
xlabel('x (1000 km)','interpreter','latex')
title('(a) Warm Pool Setup: $s^{q}, s^{\theta} (K.day^{-1})$','interpreter','latex')
ylim([0 2]);
xlim([0 40])
%  stop
end


% Plot profiles sq, so: Walker Circulation
if 1==1
subplot(iy,ix,2);
sq=          2.2*(1+0.6*cos(2*pi/LL*xx)'); % warm pool
so= 2.2*(1+0.6*cos(2*pi/LL*xx-0.1)');% warm pool walker circulation
plot(xg,oa/ta/0.1*oneday*sq,'k-'); hold on;
plot(xg,oa/ta/0.1*oneday*so,'r--'); hold on;
plot([xg(28) xg(28)],[0,2],'k--','linewidth',2);
plot(xg,oa/ta/0.1*oneday*(sq-QQ*so)/(1-QQ),'b--'); hold on;
xlabel('x (1000 km)','interpreter','latex')
title('(b) Walker Circulation Setup: $s^{q}, s^{\theta} (K.day^{-1})$','interpreter','latex')
ylim([0 2]);
xlim([0 40])
%  stop
end

% Ocean
% Plot profile eta
if 1==1
subplot(iy,ix,3)
eta = 1.5 + 0.5 *tanh(7.5*(xxo-LLo/2)); eta=eta';% eta in ocean 
%  eta=Ta/ta/0.1*oneday/Ha*eta; % optional: put dimensional
plot(xgo,eta,'k-'); hold on;
%  plot([xg(28) xg(28)],[0,0.2],'k--','linewidth',2);
xlabel('x (1000 km)','interpreter','latex')
xlim([min(xgo), max(xgo)])
%  ylim([min(eta) max(eta)]);
ylim([1 2]);
%  title('(c) Thermocline Feedback: $\eta (K.day^{-1}.m^{-1})$','interpreter','latex')
title('(c) Thermocline Feedback: $\eta(x)$','interpreter','latex')
end


% Plot gamma distribution
if 1==1
subplot(iy,ix,4)
k_sde=2;
sqref=2.2*1.6;% take at warm pool 
meanAgamma=sqref/HH; % mean value
theta_sde=(1/k_sde)*meanAgamma;
%  lambda_sde=1/(30*oneday)*ta; 
xp=(0:0.01:1)*0.4;
pass=gamma(k_sde);
Pp=1/(theta_sde(1)^k_sde)/pass*(xp.^(k_sde-1)).*exp(-xp/theta_sde(1));
plot(xp,Pp,'k-');
hold on; plot([meanAgamma meanAgamma],[0,10],'k--');
xlim([0 0.4])
ylim([0 8])
%  title('(d) Gamma Distribution: $P(a^{\prime})$','interpreter','latex');
title('(d) Gamma Distribution','interpreter','latex');
%  xlabel('$a^{\prime}$','interpreter','latex');
xlabel(' ','interpreter','latex');
%  stop;% cannot continue here
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
