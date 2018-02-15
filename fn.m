function f = fn( x, a )
% x - Variables (including both constraint and dynamic)
% a - parameters

% Create specific humidity matrix
q = zeros(1,1,2);
q(1) = x(1);
q(2) = x(2);

% Create hb and ptrop values
y = zeros(2,1);
y(1) = x(3);
y(2) = x(4);

% Run functions
F = myF( q, y, a);
G = myG( q, y, a);

% Create output vector from output of myF and myG
f = zeros(4,1);
f(1) = F(1);
f(2) = F(2);
f(3) = G(1);
f(4) = G(2);

return
