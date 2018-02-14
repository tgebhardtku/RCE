function [x,a,Del,v]=predcorr(x,a,Del,tol)

m = length(x); % number of variables
n = length(a); % number of parameters
J = zeros(m,m+n);

%Predictor
J(1:m,1:m)=fxn(x,a); 
J(1:m,m+1:m+n)=fan(x,a); 

v=null(J);

%Decide whether to use v or -v here.
if ( a * v(1) ) > 0 
   v = -v;
end

x1 = x + Del*v(1);
a1 = a + Del*v(2);

f=fn(x1,a1); 

k = 0; kmax = 7;
while ( k <= kmax ) && ( norm(f) > tol ), % stationary Newton
    %Make J have 2 components: Could move this stuff inside the loop. DONE!
    J(1:m,1:m) = fxn( x, a ); 
    J(1:m,m+1:m+n) = fan( x, a ); 
    %Use pinv to form s:
    s = -pinv(J)*f;
    x1 = x1 + s(1:m); 
    a1 = a1 + s(m+1:m+n);
    f = fn( x1, a1 ); 
    k = k + 1;
end

Del=2^((4-k)/3)*Del; % next step to be tried

% corrector succeeded, so adjust values
if ( norm(f) < tol )
    x=x1; a=a1; 
end
