%
%      Two-layer RCE model with water-vapor dependent emissivity
%
%          K. Emanuel, June 2014
%
%

function [q,y] = rce(Ts)
% Sea surface temperature, atmospheric levels, and initial humidities

% Ts     % Surface temperature (K)
Delta_T=1;    % Temperature decrease (K) from SST to 10 meters
eprecip=0.8;  % Precipitation efficiency
V=8;          % Surface wind speed (m/s)
c=0.2;        % Coefficient modifying convective MSE flux
%
%  Other parameters
%
Ck=1.5e-3;    % Surface enthalpy exchange coefficient
emmin=0.2;    % Minimum value of longwave emissivity in lower layer
%----------------------------------------------------------------------
%
%  First guess at RCE state
%
p1=750;       % First guess pressure (hPa) at first level  
pm=500;       % First guess Pressure (hPa) at mid-level    
p2=300;       % First guess Pressure (hPa) at second level 
ptrop=150;    % First guess at tropopause pressure (hPa)
rh0=0.7;      % First guess relative humidity at 10 meters
rh1=0.8;      % First guess RCE relative humidity of layer 1
rhm=0.8;      % First guess RCE relative humidity of middle level
rh2=0.8;      % First guess RCE relative humidity of layer 2
%
%   Constants
%
sig=5.67e-8;  % Stefan-Boltzmann constant
Lv=2.5e6;     % Latent heat of vaporization
cp=1005;      % Heat capacity at constant pressure
Rd=287;       % Gas constant for dry air
Rv=461;       % Gas constant for water vapor
g=9.8;        % Acceleration of gravity
pref=1000;    % Reference pressure for entropy
Tref=290;     % Reference temperature for entropy
p0=1000;      % Surface pressure (mb)
%-----------------------------------------------------------------
%
%  Normalized surface flux
%
Vnorm=1.1*Ck*V;
%
%  Find saturation moist static energy of surface
%
Tsc=Ts-273.15;
es=6.112*exp(17.67*Tsc/(243.5+Tsc));
qs=0.622.*es/(p0-es);
hs=cp*Ts+Lv*qs;
%
%  Initial guess of boundary layer moist static energy
%
hb=cp*(Ts-Delta_T)+rh0*Lv*qs;  % Boundary layer moist static energy
%
%  Iterate to find RCE boundary layer moist static energy, boundary layer
%  relative humidity, layer specific humidities, emissivities, heating
%  rates, altitudes of pressure levels, densities, dry static stabilities,
%  alpha, gamma, and convective updraft mass flux
%
hb1=hb;
hb2=hb;
delrh=1;
dtrh=0.05;
n_iteration=0;
%
while abs(delrh) > 0.001*dtrh
    n_iteration=n_iteration+1;
    if n_iteration > 1000
        disp('Warning: Maximum allowed iterations exceeded')
        disp(' ')
        break
    end    
    %    
    %  Temperatures, specific humidities, altitudes, densities, stabilities
    %
    s0=cp*log((Ts-Delta_T)/Tref)+Lv.*rh0.*qs/(Ts-Delta_T)-rh0.*qs.*Rv.*log(rh0); %PBL moist entropy
    %
    [ T1,qs1 ] = sinvert( s0,p1 );
    [ T2,qs2 ] = sinvert( s0,p2 );  % sinvert inverts saturation entropy and pressure to find temperature
    [ Tm,qsm ] = sinvert( s0,pm );
    [ Ttrop,qstrop ] = sinvert( s0,ptrop );
    %    
    z1=(hb-cp*T1-Lv*qs1)/g;
    zm=(hb-cp*Tm-Lv*qsm)/g;   % Find altitudes using conservation of moist static energy
    z2=(hb-cp*T2-Lv*qs2)/g;
    ztrop=(hb-cp*Ttrop-Lv*qstrop)/g;
    %
    q1=rh1.*qs1;              % Specific humidities
    q2=rh2.*qs2;    
    %
    S1=cp*(Tm-Ts)+g*zm;       % Dry static stabilities
    S2=cp*(Ttrop-Tm)+g*(ztrop-zm);
    %
    %  Emissivities
    %
    em1=emissivity2(p1,q1,qs1);
    em1=max(em1,emmin);  % "emissivity2" finds longwave emissivity as a function of specific humidity
    em2=emissivity2(p2,q2,qs2);
    %
    Fs=Vnorm*(hs-hb);  %  Surface turbulent enthalpy flux
    %
    % Find new estimates of effective emission temperature, tropopause
    % temperature, tropopause pressure, and layer pressures
    %
    Te=(Ts^4+Fs/sig-em1*T1^4-em2*(1-em1)*T2^4)^0.25;
    Ttrop=2^-0.25*Te;
    ptrop=pref*exp((cp*log(Ttrop/Tref)-s0)/Rd);
    pm=0.5*(p0+ptrop);
    p1=0.5*(p0+pm);
    p2=0.5*(pm+ptrop);    
    %
    % Heating rates
    %
    Q1=sig*em1*(Ts^4-2*T1^4+em2*T2^4);
    Q2=sig*em2*((1-em1)*Ts^4+em1*T1^4-2*T2^4);
    %
    alpha=Q2/(c*(Q1+Q2));
    gamma=eprecip*Q2*S1/(Q1*S2);
    delhb=hs+(Q1+Q2)/Vnorm-hb;
    hb3=hb1+dtrh*delhb;
    hb1=hb2+0.2*(hb1+hb3-2.*hb2);
    hb2=hb3;
    hb=0.5*(hb1+hb2);
    rhnew=(hb-cp*(Ts-Delta_T))/(Lv*qs);
    delrh=abs(rhnew-rh0);
    rh0=rhnew;
    %
    hm=hb-S2*gamma/c;
    qm=(hm-cp*Tm-g*zm)/Lv;
    qm=max(qm,0);
    rhm=qm/qsm;
    rh1=rhm;     %   Here we assume that the relative humidites of layers 1 and 2 are equal
    rh2=rhm;
    %
end
M=-Q1/(eprecip*S1);  % Cloud base convective updraft mass flux

q = zeros(1,1,2);
q(1) = q1;
q(2) = q2;

y = zeros(2,1);
y(1) = hb;
y(2) = ptrop;
%
%  Write output
%
% disp(['     Control Parameters']) 
% disp(['     ------------------']) %#ok<*NBRAK>
% disp(['     SST = ',num2str(Ts)])
% disp(['     V = ',num2str(V)])
% disp(['     ep = ',num2str(eprecip)])
% disp(' ')
% disp(['     Results'])
% disp(['     ------------------']) %#ok<*NBRAK>
% disp(['     Surface RH = ',num2str(rh0)])    
% disp(['     Mid-level RH = ',num2str(rhm)])  
% disp(['     Layer 1 emissivity = ',num2str(em1)])  
% disp(['     Layer 2 emissivity = ',num2str(em2)])
% disp(['     Alpha = ',num2str(alpha)])
% disp(['     Gamma = ',num2str(gamma)])
% disp(['     Layer 1 T = ',num2str(T1)])
% disp(['     Layer 2 T = ',num2str(T2)])
% disp(['     Effective emission T = ',num2str(Te)])
% disp(['     Tropopause T = ',num2str(Ttrop)])
% disp(['     Tropopause p = ',num2str(ptrop), ' mb'])
% disp(['     M= ',num2str(100*M/1.1),' cm/s'])
%