function J = fdjac(fn,x,a)
%Partial derivative with respect to x.
N = length(x);

J = eye(N);
epsilon=sqrt(eps(1));

fval = fn(x,a);
epsilon = max(1,norm(fval))*epsilon;

for j=1:N
   v = zeros(N,1);
   v(j)=epsilon;
   v = x + v;
   J(:,j) = (fn(v,a)-fval)/epsilon;
end

return
