function [ em ] = emissivity2(p, q, qs )
% Calculate infrared emissivity of thick layer as a function of water vapor
% content (in g/g), pressure 9hPa), and/or saturation specific humidity
% (g/g)
%
rh=q./qs;
em=0.5*(1+tanh((q-0.007)./0.006));
%em=q./0.01+rh;
%if p < 500
%    em=rh+0.02;
%end    
em=min(em,1);
%em=0.3;
%
end

