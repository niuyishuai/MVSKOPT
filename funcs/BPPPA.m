function xopt = BPPPA(q,p)
% q is a vector and p is a nonnegative scalar.
% solving the optimization problem
% min x'*x/2 + q'*x
% s.t. e'*x=p, x>=0.
n=length(q);
F=1:n;
xopt=zeros(n,1);
while (true)
    phi = -(p+sum(q(F)))/length(F);
    v = -(q(F)+phi);
    H = F(v<0);
    if isempty(H)
        xopt(F) = v;
        return;
    end
    F = setdiff(F,H);
end
end