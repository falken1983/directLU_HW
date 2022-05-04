function x = modbrennansc(b,a,c,f,xmin)
% Brennan-Schwartz modified Algorithm 
% Projected tridiagonal system solver with Thomas algorithm (LU factorization)
% for American Options
%
% FUNCTION : x = modbrennansc(b,a,c,f,xmin)
%
% a: Vector containing the main diagonal of the matrix
% b: Vector containing the lower diagonal of the matrix
% c: Vector containing the upperdiagonal of the matrix
% f: Vector containing the second member of the linear system
% x: Vector containing the solution of the system
% x_min : Lower obstacle in the inequality constraint

n = length(a);
% LU factorization of the matrix
alfa = zeros(1,n);
beta = zeros(1,n-1);
alfa(1) = a(1);
for i = 1:(n-1)
    beta(i) = b(i) ./ alfa(i);
    alfa(i+1) = a(i+1) - beta(i).*c(i);
end
% Initialization of solutions of the auxiliar triangular tridiagonal systems
y = zeros(1,n);
x = zeros(1,n);
% Solution of the auxiliar lower triangular tridiagonal system
y(1) = f(1);
for i = 2:n
  y(i) = f(i) - beta(i-1).*y(i-1);
end
% Solution of the auxiliar upper triangular tridiagonal system
x(n) = y(n) ./ alfa(n);
x(n) = max(xmin(n),x(n));
for i=(n-1):-1:1
  x(i) = (y(i)- c(i).*x(i+1)) ./ alfa(i);
  x(i) = max(x(i),xmin(i));
end