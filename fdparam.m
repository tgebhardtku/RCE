function fdp = fdparam(fn,x,a)
%Partial derivative with respect to parameter a.

epsilon=sqrt(eps(1));

fval = fn(x,a);
epsilon = max(1,norm(fval))*epsilon;
a1 = a + epsilon;
fdp = (fn(x,a1)-fval)/epsilon;

return
