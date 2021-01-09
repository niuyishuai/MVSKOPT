function P = genRandPortfolio(n,T)
% randomly generate portfolio datas
% n: number of assets
% T: number of periods
%% Reset random stream for reproducibility.
%rng(0,'twister');
%% Create Random Data
% Generate returns of all n assets for all T periods (between -0.1 and 0.4).
a = -0.1; b = 0.4;
D = a + (b-a)*rand(T,n);
% Generate mean returns
mu = mean(D);
% Generate covariance matrix
A = D - ones(T,1)*mu;
Sigma = (A'*A)/(T-1);
% co-skewness and co-kurtosis tensors will be created in using
P.n=n;
P.T=T;
P.A=A; % for constructing co-skewness and co-kurtosis
P.mu=mu; % mean returns
P.Sigma=Sigma; % cov
P.D=D; % raw data with all returns for all periods
end