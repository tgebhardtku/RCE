function [dely] = myG(q,y,lambda)
% this pairs with the function myF to define the DAE describing the
% two-layer model interaction for Emanuel's 2014 model.
%
%
% Inputs: q - N x N x 2 matrix with specific humidities
%         y - 2 x 1 matrix with values for hb and ptrop
%         lambda - The parameter Ts, surface temperature.
%
% Outputs: The change in y.
q = 1e-8*q;
% Constants ------------------------------------------------------------
Delta_T = 1;    % Temperature decrease (K) from SST to 10 meters

emmin=0.2;    % Minimum value of longwave emissivity in lower layer

Ck=1.5e-3;    % Surface enthalpy exchange coefficient
V = 8;
Vmin = 2;

sig=5.67e-8;  % Stefan-Boltzmann constant
Lv=2.5e6;     % Latent heat of vaporization
cp=1005;      % Heat capacity at constant pressure
Rd=287;       % Gas constant for dry air
Rv=461;       % Gas constant for water vapor
pref=1000;    % Reference pressure for entropy
Tref=290;     % Reference temperature for entropy
p0=1000;      % Surface pressure (mb)

% Fixed variables ------------------------------------------------------
hb = y(1);
ptrop = y(2);

Ts = lambda(1);

% Calculations for y variables -------------------------------------
% Surface wind normalization
Vmin = Vmin*1.1*Ck;
Vnorm = V*1.1*Ck;

%  Find saturation moist static energy of surface
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

% "emissivity" finds longwave emissivity as a function of specific humidity
% NOTE: In this case, we will treat the true emmissivity as the average over
% space.
em1=emissivity2(p1,q(:,:,1),qs1);
em1=max(em1,emmin);  
em1=mean(mean(em1));
em2=emissivity2(p2,q(:,:,2),qs2);
em2=min(em2,1);
em2=mean(mean(em2));

%  Radiative heating rates in each layer
Q1=sig*em1.*(Ts.^4-2.*T1.^4+em2.*T2^4);
Q2=sig*em2.*((1-em1).*Ts^4+em1.*T1^4-2.*T2^4);

% Surface enthalpy flux 
Vs=Vnorm;
Vs=max(Vs,Vmin);
Fs=Vs.*(hs-hb);

% Update calculations --------------------------------------------------

hb_new = hs+(Q1+Q2)/Vnorm;

Te = (Ts^4+Fs/sig-em1*T1^4-em2*(1-em1)*T2^4)^0.25;
Ttrop = 2^-0.25*Te;
ptrop_new = pref*exp((cp*log(Ttrop/Tref)-s0)/Rd);

dely(1) = hb_new - hb;
dely(2) = ptrop_new - ptrop;