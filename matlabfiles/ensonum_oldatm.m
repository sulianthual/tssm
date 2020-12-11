 % by Sulian Thual
%
% ENSO model with WWB, final and clean version
% Function, Solves generic atmospheric model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ka=ensonum_oldatm(schema,nxa,dx,ddk,cck,HK)
%
% Generic model: 
% ddk*Ka + cck* dKa/dx = HK-<HK>
% with periodic boundary conditions
%
% Given HK, ddk,cck , returns K
%
% nx=size of grid (entire equatorial belt)
% < >= zonal average (entire equatorial belt)
% dx=grid step
%
% Solving is performed using AX=B, X=(A^-1)B:
% which gives K=X, where A=ddk(..)+cck*d(..)/dx and B=HK-<HK>
% The dissipation ddk cannot be zero or matrix A is singular (non-invertible)
% 
%%%%%%%%%%%%%%%%%%%%%
%
%  % Method centered difference dK/dx=(Ki+1-Ki-1)/2dx: AVOID THIS SCHEME (unusable directly due to periodic boundaries)
%  if schema==1
%  % Compute matrix A
%  AA=zeros(nxa,nxa);
%  for ix=1:nxa-1; AA(ix,ix+1)= cck/dx/2; end; AA(nxa,1)=cck/dx/2; % A(i,i+1) 
%  for ix=1:nxa;   AA(ix,ix)= ddk ; end
%  for ix=2:nxa;   AA(ix,ix-1)= -cck/dx/2; end; AA(1,nxa)= -cck/dx/2;  % A(i,i-1) 
%  % Compute vector B
%  BB=zeros(nxa,1);
%  for ix=1:nxa; BB(ix,1)= HK(ix,1); end; 
%  BB(1:nxa,1)=BB(1:nxa,1)-mean(BB(1:nxa,1));% remove zonal average
%  % Compute X (matrix inversion) 
%  XX=AA\BB;
%  Ka=XX(1:nxa); 
%  end

% Method forward Euler dK/dx=(Ki+1-Ki)/dx: AVOID (slightly less accurate)
if schema==2
% Compute matrix A
AA=zeros(nxa,nxa);
for ix=1:nxa;   AA(ix,ix)= ddk +cck/dx ; end
for ix=2:nxa;   AA(ix,ix-1)= -cck/dx; end; AA(1,nxa)= -cck/dx;  % A(i,i-1) 
% Compute vector B
BB=zeros(nxa,1);
for ix=1:nxa; BB(ix,1)= HK(ix,1); end; 
BB(1:nxa,1)=BB(1:nxa,1)-mean(BB(1:nxa,1));% remove zonal average
% Compute X (matrix inversion) 
XX=AA\BB;
Ka=XX(1:nxa); 
end

% Method backward Euler dK/dx=(Ki-Ki-1)/dx: AVOID (slightly less accurate)
if schema==3
% Compute matrix A
AA=zeros(nxa,nxa);
for ix=1:nxa-1; AA(ix,ix+1)= cck/dx; end; AA(nxa,1)=cck/dx; % A(i,i+1) 
for ix=1:nxa;   AA(ix,ix)= ddk -cck/dx; end
% Compute vector B
BB=zeros(nxa,1);
for ix=1:nxa; BB(ix,1)= HK(ix,1); end; 
BB(1:nxa,1)=BB(1:nxa,1)-mean(BB(1:nxa,1));% remove zonal average
% Compute X (matrix inversion) 
XX=AA\BB;
Ka=XX(1:nxa); 
end

%  


%%%%%%%%%%%%%%%%%%%%%
%
