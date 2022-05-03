% CN Survey
% Structural Data
% a = 0.03; sigma = 0.0075; K = .975;
% K = 0.95;
% vertices = [.5 1:10];
% zeroes = 0.01*[.03 .03 .12 .27 .47 .68 .91 1.13 1.35 1.54 1.73];
% T = 2; Tau =5;
% T = 3; Tau =5;

% Structural Data Thesis
a = 0.1; sigma = 0.01; K = .85;
vertices = [.083 .25 1:10];
zeroes = 0.01*[3.46 3.54 4.02 4.51 4.79 4.98 5.13 5.24 5.35 5.44 5.51 5.56];
T = 3; Tau =5;
% % T = 2; Tau =5;

opt = 'us';   % American Options
interp ='hw'; % TSIR Estimation w/ Minimal Consistent Family 

steps =[.5 .2 .1 .05 .02 .01 .005 .002 .001 .0005]; 
% steps =[.5 .2 .1 .05 .02 .01 .005 .002 .001]; 
L = 10000;

N = size(steps,2);

% Preallocating Prices
p = zeros(N, 4);
elapsed_time = zeros(N, 4);
% errp = zeros(N-1,4);

% pSOR Solver Out Data
miniter = zeros(N,2);
medianiter = zeros(N,2);
maxiter = zeros(N,2);
om = zeros(N,2);

% CN (pSOR Algo)%
disp('SQCN w/ pSOR Algo')
for i = 1:N
    disp(steps(i))
    t1 = tic; [price, size_lattice, iter, omega] = ...
        CNMixedBondOptionOPT2(a, sigma, vertices, zeroes, steps(i), T, Tau, K, interp, opt, 'efd', 1.); ...
        elapsed_time(i,1)=toc(t1); p(i,1) = L*price; 
    maxiter(i,1) = max(iter); miniter(i,1) = min(iter); medianiter(i,1) = median(iter);
    om(i,1) = omega;
end

% CN++ (pSOR Algo)% 
disp('LCN w/ pSOR Algo')
for i = 1:N
    disp(steps(i))
    t1 = tic; [price, size_lattice, iter, omega] = ...
        CNMixedBondOptionOPT2(a, sigma, vertices, zeroes, steps(i), T, Tau, K, interp, opt, 'cn', sqrt(2.)); ...
        elapsed_time(i,2)=toc(t1); p(i,2) = L*price; 
    maxiter(i,2) = max(iter); miniter(i,2) = min(iter); medianiter(i,2) = median(iter);
    om(i,2) = omega;
end

% SQCN, B-S Algo (CN)
disp('SQCN w/ B-S Algo')
for i = 1:N
    disp(steps(i))
    t1 = tic; price = ...
        CNMixedBondOptionOPT3(a, sigma, vertices, zeroes, steps(i), T, Tau, K, interp, opt, 'efd', 1.); ...
        elapsed_time(i,3) = toc(t1); p(i,3) = L*price;
end

% LCN, B-S Algo (CN++)
disp('LCN w/ B-S Algo')
for i = 1:N
    disp(steps(i))
    t1 = tic; price = ...
        CNMixedBondOptionOPT3(a, sigma, vertices, zeroes, steps(i), T, Tau, K, interp, opt, 'cn', sqrt(2.)); ...
        elapsed_time(i,4) = toc(t1); p(i,4) = L*price;
end

errp = abs(diff(log(p),1,1));
% EFD
% for i = 1:N
%     steps(i)
%     t1 = tic; price = ...
%         UEFDMBondOptionOPT2(a, sigma, vertices, zeroes, steps(i), T, Tau, K, interp, opt, 1./3.); ...
%         toc(t1); p(i) = L*price;    
% end