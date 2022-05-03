function [R,pesos,AD,Omega] = trinHW(T,N,a,sigma,P)
%  trinHW: Genera un árbol de Hull y White
%  Sintaxis: 
%
%            arbol = trinHW(T,N,a,sigma,P)
%
%        T     : Tiempo de vencimiento de la opción
%        N     : Número de pasos
%        a     : Tasa de reversión
%        sigma : Volatilidad
%        P     : Precios de los bonos en los tiempos del árbol
%
%        R     : Árbol de tipos 
%        pesos : Pesos de cada trayectoria
%        AD    : Precio Arrow-Debreu
%        Omega : Auxiliar
%
%        Ejemplo:   [R,pesos,AD,Omega] = arbolHW(1,2,0.1,0.01,[0.95123 0.88692])
%
global L zeroindex dt nfilas AnalyticOmega
%
%    Parámetros del árbol
%
dt = T/N;                % Intervalo de tiempo en árbol
%
M = -a*dt;
V = sigma^2*dt;
H = sqrt(3*V);
L = fix(0.1835/abs(M));
%
%
zeroindex = L+3;
nfilas = 2*(L+2)+1;
pesos = zeros(nfilas,3);
St = H * ((L+2):-1:-(L+2))';
%
%    Pesos: Configuración central
%
indices = (-L:L)';
auxindices = zeroindex-indices;
indices2 = indices.*indices;
pesos(auxindices,1) = 1.0/6.0 + 0.5*(indices2*M^2 + indices*M); % up
pesos(auxindices,2) = 2.0/3.0 - indices2*M^2;                   % middle
pesos(auxindices,3) = 1.0/6.0 + 0.5*(indices2*M^2 - indices*M); % down
%
%    Pesos: Configuración descendente
%
indices = L+1;
indices2 = indices.*indices;
pesos(2,2) =  7.0/6.0 + 0.5*(indices2*M^2 + 3*indices*M);  % middle
pesos(2,3) = -1.0/3.0 - indices2*M^2-2*indices*M;         % down
pesos(2,1) =  1.0/6.0 + 0.5*(indices2*M^2 + indices*M);    % downdown
%
%    Pesos: Configuración ascendente
%
indices = -(L+1);
indices2 = indices.*indices;
pesos(nfilas-1,3) =  1.0/6.0 + 0.5*(indices2*M^2 - indices*M);    % upup
pesos(nfilas-1,1) = -1.0/3.0 - indices2*M^2+2*indices*M;          % up
pesos(nfilas-1,2) =  7.0/6.0 + 0.5*(indices2*M^2 - 3*indices*M);  % middle
%
%    Precio Arrow-Debreu
%
AD = zeros(nfilas,N);
AD(zeroindex,1) = 1;
Disct = zeros(nfilas,N);
R = zeros(nfilas,N);
R(zeroindex,1) = -log(P(1))/dt;   % Tipo de interés para el primer intervalo
Omega(1) = -log(P(1))/dt;
% Omega(1)=AnalyticOmega(1);
for n = 2: N
   auxn = min(n-1,L+1);
   indices = -auxn:auxn;
   auxindices = zeroindex-indices;
   Disct(:,n-1) = exp(-R(:,n-1)*dt); 
   AD(auxindices,n) = pesos(auxindices-1,3).*Disct(auxindices-1,n-1).*AD(auxindices-1,n-1) +... %up
      pesos(auxindices,2).*Disct(auxindices,n-1).*AD(auxindices,n-1) +... % middle
      pesos(auxindices+1,1).*Disct(auxindices+1,n-1).*AD(auxindices+1,n-1);  % down
   %
   %    Configuración ascendente
   %
   AD(nfilas-3,n) = AD(nfilas-3,n) + pesos(nfilas-1,3) * Disct(nfilas-1,n-1)*AD(nfilas-1,n-1);
   %
   %    Configuración descendente
   %
   AD(4,n) = AD(4,n) + pesos(2,1) * Disct(2,n-1)*AD(2,n-1);
   %
   %   Cálculo de omega y del tipo de interés
   %
   Omega(n) = log(exp(-dt*St)'*AD(:,n)/P(n))/dt;
% Omega(n)=AnalyticOmega(n);
R(auxindices,n) =  St(auxindices) + Omega(n);
end
