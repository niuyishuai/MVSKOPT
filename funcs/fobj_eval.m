function [fval,dfval] = fobj_eval(x,c,P,opt)
    % evaluate MVSK polynomial objective function and gradient at given point x
    % f(x) = - c1*(mu'*x) + c2*(x'*Sigma*x) 
    %        - c3*x'*coskewness*kron(x,x) + c4*(x'*ckkurtosis*kron(kron(x,x),x)
    % df(x) = - c1*mu + 2*c2*Sigma*x 
    %        - 3*c3*coskewness*kron(x,x) + 4*c4*ckkurtosis*kron(kron(x,x),x)
    xx=kron(x,x);
    xxx=kron(xx,x);
    A = P.Sigma*x;
    B = P.coskewness*xx;
    C = P.cokurtosis*xxx;
    if opt==1 % evaluate function only
        fval = -c(1)*P.mu'*x+c(2)*x'*A-c(3)*x'*B+c(4)*x'*C;
    else % evaluate function and gradient
        fval = -c(1)*P.mu'*x+c(2)*x'*A-c(3)*x'*B+c(4)*x'*C;
        dfval = -c(1)*P.mu + 2*c(2)*A - 3*c(3)*B + 4* c(4)*C;
    end
end