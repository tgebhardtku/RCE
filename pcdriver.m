tol=1.E-10;
max_it = 80;
x = zeros(4,1);
xarr = zeros(4,max_it+1);
SSTarr = zeros(1,max_it+1);
jcond = zeros(1,max_it); % Records the conditioning of the jacobian
SST = 300; %Surface temperature (K)
SSTarr(1) = SST;
% Converge to equilibrium
[q,y] = rce(SST);
x(1) = 1e8*q(1); % First level specific humidity (variable)
x(2) = 1e8*q(2); % Second level specific humidity (variable)
x(3) = y(1); % Boundary layer moist static energy (constraint)
x(4) = y(2); % Pressure at tropopause (constraint)
xarr(:,1) = x;
% Use predcorr to continue equilibrium
Del=0.1;
for k=1:max_it
[x,SST,Del,v]=predcorr(x,SST,Del,tol);
jcond(k) = cond(fxn(x,SST));
xarr(:,k+1)=x;
SSTarr(k+1)=SST;
disp(k);
end

% figure(1);
% plot(xarr,aarr,'o-')
%hold on
%pause(.01);
%figure(2)
%plot(k,v(1),'+',k,v(2),'x');
%hold on
