function [ T,qs] = sinvert( s0,p )
%  This function accepts saturation entropy, s0, and pressure 
% (hPa), and returns temperature (K), and saturation
% specific humidity (g/g)
%  Note:  This is not a very efficient algorithm
%
g=9.8;
cp=1005;
Rd=287;
Lv=2.5e6;
del=10;
pref=1000;
Tref=290;
% tolerance=.002;
tolerance=2e-8;
%
%   First guess
%
Tg=280;
while del > tolerance
    %
    Tgc=Tg-273.15;
    esg=6.112*exp(17.67*Tgc/(243.5+Tgc));
    qsg=0.622*esg/(p-esg);
    T=Tref*exp((s0+Rd*log(p/pref)-Lv*qsg/Tg)/cp);
    del=abs(T-Tg);
    Tg=0.65*Tg+0.35*T;
end
T=Tg;
Tgc=Tg-273.15;
esg=6.112*exp(17.67*Tgc/(243.5+Tgc));
qs=0.622*esg/(p-esg);
end

