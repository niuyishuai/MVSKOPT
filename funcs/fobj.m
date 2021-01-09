function [f]=fobj(x,c,m1,m2,m3,m4)
% c is investor's preference level with sum(c)=1, c>=0
% x decision variable
f = -c(1)*m1+c(2)*m2-c(3)*m3+c(4)*m4;
end