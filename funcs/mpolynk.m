function p=mpolynk(n,k)
% generate a polynomial of n variables and k monomials
p=mpoly(n);
p.k=k;
p.pow=zeros(k,n);
p.coef=zeros(k,1);
end