function P = genRandPortfolio_mo(n,T)
% randomly generate portfolio datas
% n: number of assets
% T: number of periods
%% Reset random stream for reproducibility.
%rng(0,'twister');
%% Create Random Data
% Generate returns of all n assets for all T periods (between -0.1 and 0.4).
a = -0.1; b = 0.4;
D = a + (b-a)*rand(T,n);
% Generate mean, variance, co-skewness and co-kurtosis trnsors
mu = mean(D);
coskewness=zeros(n,n^2);
cokurtosis=zeros(n,n^3);
s_demeaned=D(:,:)-kron(mu,ones(T,1)); %creates a matrix T * N with generic element [a_ij-mu_j] with average mu_j of the j_th asset

varcov=1/T*(s_demeaned)'*(s_demeaned);
for i=1:T
    coskewness=coskewness+1/T*(kron(s_demeaned(i,:)'*s_demeaned(i,:),s_demeaned(i,:)));
    cokurtosis=cokurtosis+1/T*(kron(s_demeaned(i,:)',kron(s_demeaned(i,:),kron(s_demeaned(i,:),s_demeaned(i,:)))));
end
P.n=n;
P.T=T;
P.mu=mu(:); % mean
P.Sigma=varcov; %variance
P.coskewness=coskewness; %coskewness
P.cokurtosis=cokurtosis; %cokurtosis
P.D=D; % raw data with all returns for all periods
%% Generate MVSK model constraint
P.Cons.Aeq = ones(1,n);
P.Cons.beq = 1;
P.Cons.lb = zeros(n,1);
P.Cons.ub = ones(n,1);
end