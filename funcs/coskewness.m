function f = coskewness(i,j,k,P)
% computing co-skewness S(i,j,k)
f=sum(P.A(:,i).*P.A(:,j).*P.A(:,k))/P.T;
end