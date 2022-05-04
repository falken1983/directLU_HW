function x = tridiagsys(b,a,c,f)
% Tridiagonal system solver with Thomas algorithm (LU factorization)
%
% FUNCTION : x = tridiagsys(b,a,c,f)
%
% a: Vector containing the main diagonal of the matrix
% b: Vector containing the lower diagonal of the matrix
% c: Vector containing the upperdiagonal of the matrix
% f: Vector containing the second member of the linear system
% x: Vector containing the solution of the system

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
for i=(n-1):-1:1
  x(i) = (y(i)- c(i).*x(i+1)) ./ alfa(i);
end