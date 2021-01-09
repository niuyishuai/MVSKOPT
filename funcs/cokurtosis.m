function f = cokurtosis(i,j,k,l,P)
% computing co-skewness K(i,j,k,l)
f=sum(P.A(:,i).*P.A(:,j).*P.A(:,k).*P.A(:,l))/P.T;
end