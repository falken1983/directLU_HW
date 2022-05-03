% Structural Data
a = 0.03; sigma = 0.0075; K = .975;
vertices = [.5 1:10];
zeroes = 0.01*[.03 .03 .12 .27 .47 .68 .91 1.13 1.35 1.54 1.73];
T = 3; Tau =5;
opt = 'us';   % American Options
interp ='hw'; % TSIR Estimation w/ Minimal Consistent Family 

% steps =[.5 .2 .1 .05 0.02 0.01 0.005 .002 .001 .0005];

steps =[.5 .2 .1 .05 0.02 0.01 0.005 .002 .001];

% steps =[.5 .2 .1 .05 0.02 0.01 0.005];
% steps = [.005 .002 .001];
% steps =[.5 .2 .1 .05 ]; 
% steps =.05; 
% steps = .001;

L = 10000;

% Preallocating Prices
N = size(steps,2);
p = zeros(1, N);
elapsed_time = zeros(1,N);

miniter = zeros(N,1);
medianiter = zeros(N,1);
maxiter = zeros(N,1);
om = zeros(N,1);
% CN
% for i = 1:N
%     steps(i)
%     t1 = tic; [price, size_lattice, iter, omega] = ...
%         CNMixedBondOptionOPT2(a, sigma, vertices, zeroes, steps(i), T, Tau, K, interp, opt, 'efd', 1.); ...
%         toc(t1); p(i) = L*price; 
%     maxiter(i) = max(iter); miniter(i) = min(iter); medianiter(i) = median(iter);
%     om(i) = omega;
% end
for i = 1:N
    steps(i)
    t1 = tic; [price, size_lattice, iter, omega] = ...
        CNMixedBondOptionOPT2(a, sigma, vertices, zeroes, steps(i), T, Tau, K, interp, opt, 'cn', sqrt(2.)); ...
        elapsed_time(i)=toc(t1); p(i) = L*price; 
    maxiter(i) = max(iter); miniter(i) = min(iter); medianiter(i) = median(iter);
    om(i) = omega;
end