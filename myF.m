function [fq] = myF(q, y, lambda)
% This is the forward operator for the two-level model from Emanuel.
% Lambda is Ts, the surface temperature. Fixed are all of the variables
% from the rce calculation that are fixed in this one.
%
%
% Inputs: q - N x N x 2 matrix with specific humidities
%         y - 2 x 1 matrix with values for hb and ptrop
%         lambda - The parameter Ts, surface temperature.
%
% Outputs: The change in q.
fq = zeros(size(q));
q = 1e-8*q;
% Constants ------------------------------------------------------------
Delta_T = 1;    % Temperature decrease (K) from SST to 10 meters
eprecip = 0.7;  % Precipitation efficiency
c = 0.5;        % Coefficient modifying convective MSE flux

emmin=0.2;    % Minimum value of longwave emissivity in lower layer

Ck=1.5e-3;    % Surface enthalpy exchange coefficient
V = 8;
Vmin = 2;

sig=5.67e-8;  % Stefan-Boltzmann constant
Lv=2.5e6;     % Latent heat of vaporization
cp=1005;      % Heat capacity at constant pressure
Rv=461;       % Gas constant for water vapor
g=9.8;        % Acceleration of gravity
Tref=290;     % Reference temperature for entropy
p0=1000;      % Surface pressure (mb)

% Fixed variables ------------------------------------------------------
hb = y(1);
ptrop = y(2);

% Calculations for fixed variables -------------------------------------
% Surface wind normalization
Vmin = Vmin*1.1*Ck;
Vnorm = V*1.1*Ck;

%  Find saturation moist static energy of surface
Ts = lambda(1);
Tsc=Ts-273.15;
es=6.112*exp(17.67*Tsc/(243.5+Tsc));
qs=0.622.*es/(p0-es);
hs=cp*Ts+Lv*qs;

% Find pressure levels
pm = 0.5*(ptrop + p0);
p1 = 0.5*(pm + p0);
p2 = 0.5*(pm + ptrop);

% PBL moist entropy calculation
rh0 = (hb-cp*(Ts-Delta_T))/(Lv*qs);
s0=cp*log((Ts-Delta_T)/Tref)+Lv.*rh0.*qs/(Ts-Delta_T)-rh0.*qs.*Rv.*log(rh0);

% Invert saturation entropy and pressure to find temperature
[ T1,qs1 ] = sinvert( s0,p1 );
[ T2,qs2 ] = sinvert( s0,p2 );
[ Tm,qsm ] = sinvert( s0,pm );

% Find altitudes using conservation of moist static energy
z2=(hb-cp*T2-Lv*qs2)/g;
zm=(hb-cp*Tm-Lv*qsm)/g;

% Dry static stabilities
S1=cp*(Tm-Ts)+g*zm;
S2=cp*(T2-Tm)+g*(z2-zm);

%  Layer capacities

c1=100*(p0-pm)*Lv/g;
c2=100*(pm-ptrop)*Lv/g;
% Update calculations --------------------------------------------------

% Update mid-level moist static energy
hm=cp*Tm+g*zm+0.5*Lv*(q(:,:,1)+q(:,:,2));

% "emissivity" finds longwave emissivity as a function of specific humidity
em1=emissivity2(p1,q(:,:,1),qs1);
em1=max(em1,emmin);  
em2=emissivity2(p2,q(:,:,2),qs2);
em2=min(em2,1);

%  Radiative heating rates in each layer
Q1=sig*em1.*(Ts.^4-2.*T1.^4+em2.*T2^4);
Q2=sig*em2.*((1-em1).*Ts^4+em1.*T1^4-2.*T2^4);

% Surface enthalpy flux 
Vs=Vnorm;
Vs=max(Vs,Vmin);
Fs=Vs.*(hs-hb);

alpha = Q2/(c*(Q1+Q2));
gamma = eprecip*Q2*S1/(Q1*S2);

delhinv=1./max(hb-hm,1e3);

% First layer vertical velocity
w1=(alpha.*eprecip.*Fs.*delhinv+Q1./S1)./(1-eprecip); 
w1=w1-mean(mean(w1));  % Force area mean w1 to be zero

% Convective updraft mass flux
M=w1+alpha.*Fs.*delhinv;
M=max(M,0);

% Second layer vertical velocity
w2=gamma.*M+Q2./S2;
w2 = w2 - mean(mean(w2)); % Force area mean w2 to be zero

% Generate output ------------------------------------------------------
fq(:,:,1)=((Fs+(w1-c*M).*(hb-hm))+Q1)./c1;
fq(:,:,2)=((c*M-w2).*(hb-hm)+Q2)./c2;
fq = 1e8*fq;