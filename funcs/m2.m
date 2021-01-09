function f=m2(x,Sigma)
% variance
% f=m2(x,Sigma)
% where Sigma is covariance matrix
f=transpose(x)*Sigma*x;
end