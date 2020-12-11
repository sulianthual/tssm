function Ke = ensonum_wavefourier(K,Ek,Ek0,Ck,Pk,nx,dx,dt)
% Solves evolution of a long-wave using analytic solution in zonal Fourier space
% Sulian Thual CIMS fev 2013
%
% solves: 
% dtK+Ek*K+Ck*dxK=Pk with periodic boundary condition
% which gives in zonal Fourier space (dx=i*k): 
% dtK+ik*Ck*K=Pk
% K(k,t+dt)=K(k,t)exp(-i*k*Ck*dt)+Pk/(i*k) FOR k ne 0
% K(0,t+dt)=K(0,t) + Pk*dt
% 
% Inputs :
% - K(k,t), in zonal Fourier space  at time t 
% - speed Ck 
% - forcing Pk(k,t) in zonal Fourier space at time t
% - k: zonal wavenumber (give row of values) (in Fourier it can also be 2*pi*k)
% - dt: timestep
% - Ek: dissipation coefficient
% - Ek0: potentially different dissipation coefficient for the k=0 component
%
% Outputs:
% K(k,t+dt)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
k=2*pi/(nx*dx)*(-nx/2:nx/2-1);% pulsation wavenumber
ii=complex(0,1); % the complex number
%
% Go to Fourier
K=fftshift(fft(K));
Pk=fftshift(fft(Pk));

% Solve Evolution exactly
Ke=K.*exp(-Ek*dt-ii*k'*Ck*dt) +(Pk./(Ek+ii*k'*Ck)).*(1-exp(-Ek*dt-ii*k'*Ck*dt));


% Solve the k=0 component differently
[i0,term]=ensospe_searchclosest(k,0); %k(i0)
if Ek0==0;% Ek0 is the damping for k=0 (potentially different from Ek)
Ke(i0)=K(i0) + Pk(i0)*dt;
else
Ke(i0)=K(i0)*exp(-Ek0*dt) +Pk(i0)/Ek0*(1-exp(-Ek0*dt));
end

% Go back to physical
Ke=ifft(ifftshift(Ke));
Ke=real(Ke); % remove potential round-off imaginary part
