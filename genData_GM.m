%% EECS 545 project
% generates clusters with Gaussian distribution plus noise for 
% testing clustering algorithms
% Charles Lu

%% Parameters
C = 4;      % number of clusters
N_c = 100;  % number of points per cluster
N_n = 100;  % number of noise points not belonging to any cluster
D = 50;      % number of dimensions

%% Process parameters
N = C*N_c + N_n;

%% Create points
X = zeros(N,D);
Y = zeros(N,1);
for c = 1:C
    mu = 10*(rand(1,D)-0.5);
    M = rand(D,C);
    sigma = (M*M'+diag(rand(D,1)));
    sigma = sigma/max(max(sigma));
    X(1+(c-1)*N_c:c*N_c,:) = mvnrnd(mu,sigma,N_c);
    Y(1+(c-1)*N_c:c*N_c) = c;
end
X(c*N_c+1:end,:) = 20*(rand(N_n,D)-0.5);
Y(c*N_c+1:end) = -1;