function P = genPortfolio(n,T,D)
% Generate mean returns
mu = mean(D);
% Generate covariance matrix
A = D - ones(T,1)*mu;
Sigma = (A'*A)/T;
% Generate co-skewness tensor
S=zeros(n,n,n);
for i=1:n
    for j=i:n
        for k=j:n
            v=sum(A(:,i).*A(:,j).*A(:,k))/T;
            idx=perms([i,j,k]); % get all permutation of [i,j,k]
            for u=idx'
                S(u(1),u(2),u(3))=v;
            end
        end
    end
end
% Generate co-kurtosis tensor
K=zeros(n,n,n,n);
for i=1:n
    for j=i:n
        for k=j:n
            for l=k:n
                v=sum(A(:,i).*A(:,j).*A(:,k).*A(:,l))/T;
                idx=perms([i,j,k,l]); % get all permutation of [i,j,k,l]
                for u=idx'
                    K(u(1),u(2),u(3),u(4))=v;
                end
            end
        end
    end
end
P.mu=mu; % mean returns
P.Sigma=Sigma; % cov
P.S=S; % coskewness
P.K=K; % cokurtosis
P.D=D; % raw data with all returns for all periods
end