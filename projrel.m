function [x,iter] = projrel(b,a,c,f,x_min,omega,max_iter,error_max,xin)
% Projected relaxation for tridiagonal systems
%
% Solution of the complementarity problem
%
% Ax <=f , x >=x_min, (Ax-f).(x-x_min) = 0
%
% FUNCTION : [x,iter] = projrel(b,a,c,f,x_min,omega,max_iter,error_max)
%
%Data parameters
%
% b a c: Vectors containing the tridiagonal matrix
% a: Main diagonal of the matrix
% b: Lower diagonal of the matrix
% c: Upper diagonal of the matrix
% f: Second member of the system
% x_min : Lower obstacle in the inequality constraint
% omega : Relaxation parameter
% max_iter : Maximum number of iterations 
% error_max : Tolerance for the stopping test in relaxation algorithm
%xin :initial iteration vector
%
% Output
%
% x: Solution
% iter : Number of iterations

n = length(a);
x = xin;
dif = zeros(1,n);
error = error_max + 1;
iter = 0;

while ((error > error_max) & (iter < max_iter))
  aux = x;
  x(1) = x(1) + omega.*(f(1) - c(1).*x(2)-a(1).*x(1))./a(1);
  x(1) = max(x(1),x_min(1)); % Blocking this line for solving Ax=f
  for i=2:n-1
    x(i) = x(i) + omega.*(f(i) - x(i-1).*b(i-1) - x(i+1).*c(i)-a(i).*x(i))./a(i);
    x(i) = max(x(i),x_min(i)); % Blocking this line for solving Ax=f
  end;
  x(n) = x(n) + omega.*(f(n) - x(n-1).*b(n-1)-a(n).*x(n))./a(n);
  x(n) = max(x(n),x_min(n)); % Blocking this line for solving Ax=f
  dif = x - aux;
  error = (sum(dif.^2))./sum(x.^2);
  iter = iter + 1;  
end;  
