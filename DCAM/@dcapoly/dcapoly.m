classdef dcapoly < matlab.mixin.Copyable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dca class for polynomial optimization (without Yalmip)
    % dc programming solver
    %
    % Author: yi-shuai niu
    % 2019-4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        dcp % dc program
        x0 % initial point
    end
    properties(GetAccess = public,SetAccess = private)         % read only
        fopt = inf % objective value
        xopt = [] % optimal solution
        iter = 0 % iterations of dca
    end
    properties
        tolf = 1e-6 % tolerence for objective function
        tolx = 1e-6 % tolerence for iterative point
        maxiter = 10000 % max iterations for dca
        plot = 0  %1: draw iterations of dca, 0: otherwise
        verbose = 1  %1: display iterations of dca, 0: otherwise
        convexsolver = 'quadprog' % convex subproblem solver, default quadprog
        convexsolver_option = optimset( 'Display','off'); %
        linesearch = false % use line search for acceleration
        plotlinetype='b-s'; % if plot = 1, this option set the plot line type
        plotinnewfig=true; % if true, we flot in a new fig, otherwise, we plot in fig 1
        rho;
        P;
        c;
        stepsize=inf; % initial step size for armijo rule
        qpinitialization=true;
    end
    
    methods
        function obj = dcapoly(dcp,x0)
            % dca constructor
            % obj = dca(dcp,x0)
            % where dcp is a dc program object
            % x0 is an initial point. It will be a random point if x0 is not given.
            if nargin==1
                obj.dcp = dcp;
                x0 = rand(size(dcp.X));
                obj.x0 = x0(:);
            elseif nargin==2
                obj.dcp = dcp;
                obj.x0 = x0(:);
            else
                error('wrong input arguments.');
            end
        end
        function xopt = get.xopt(obj)
            % get optimal solution
            xopt = obj.xopt;
        end
        function fopt = get.fopt(obj)
            % get objective value
            fopt = obj.fopt;
        end
        function set.tolf(obj,val)
            % set tolerence of objective function
            obj.tolf = val;
        end
        function set.tolx(obj,val)
            % set tolerence of iterative point
            obj.tolx = val;
        end
        function set.maxiter(obj,val)
            % set max iterations of dca
            obj.maxiter = val;
        end
        function set.plot(obj,val)
            % set plotting option of dca
            obj.plot = val;
        end
        function set.verbose(obj,val)
            % set verbose option of dca
            obj.verbose = val;
        end
        function set.convexsolver(obj,val)
            % set convex subproblem solver (used for yalmip)
            obj.convexsolver = val;
        end
        function set.linesearch(obj,yn)
            % set line search
            obj.linesearch = yn;
        end
        function set.rho(obj,val)
            % set parameter rho for universal decomposition
            obj.rho = val;
        end
        % dca algorithm for solving dcp with starting point x0
        function status = optimize(obj)
            % dca optimizer
            % status.flag : 0 dca converges with tolx or tolf
            %               1 maxiter exceed
            %               2 problem infeasible or unbounded
            % status.info : solution informations.
            % status.iter : number of iterations of dca.
            % status.time : cpu time for dca (sec.)
            % status.avgt : average time for each iteration (sec.)
            
            obj.iter = 0;
            xk = obj.x0;
            % plotting if actived
            if (obj.plot==1)
                if (obj.plotinnewfig)
                    figure
                else
                    figure(1);
                end
            end
            % display of actived
            if (obj.verbose == 1)
                fprintf('------------------------------------------------------------\n');
                fprintf('DCA version 1.0 beta \nLocal solver: %s\n',obj.convexsolver);
                if obj.linesearch
                    fprintf('* activate Armijo linesearch acceleration\n');
                end
                fprintf('------------------------------------------------------------\n');
                fprintf('Iterations | Objective values |   Delta x   |   Delta f \n');
            end
            cputime=tic;
            
            [fk,dfk]=obj.dcp.F.f(xk,0); % compute f(xk) and df(xk)
            while obj.iter < obj.maxiter
                obj.iter = obj.iter+1;
                if obj.qpinitialization && obj.iter==1
                    xk1 = quadprog(2*obj.c(2)*(obj.P.Sigma), -obj.c(1)*obj.P.mu,[],[],obj.P.Cons.Aeq,obj.P.Cons.beq,obj.P.Cons.lb,obj.P.Cons.ub,xk,obj.convexsolver_option);
                    [fk1,dfk1] = obj.dcp.F.f(xk1,0);
                else
                    yk = obj.rho*xk - dfk; % compte dh(xk)
                    xk1=BPPPA(-yk/obj.rho,1);
                    [fk1,dfk1] = obj.dcp.F.f(xk1,0);
                    % accelerate with line search
                    if obj.linesearch == true
                        [feasxk,activesetxk] = checkfeas(xk);
                        [feasxk1,activesetxk1] = checkfeas(xk1);
                        d = xk1 - xk;
                        if (feasxk && feasxk1 && min(activesetxk - activesetxk1)>=0 && d'*dfk1 < 0)
                            [xacc,facc,obj.stepsize]=armijo_adaptive(obj,fk1,xk1,d,obj.stepsize); % increase many times
                            if obj.verbose == 1
                                fprintf('accelerated: reduced %17.3e  moved %17.3e \n',facc-fk1,norm(xacc-xk1));
                            end
                            fk1 = facc;
                            xk1 = xacc;
                        end
                    end
                end
                % compute errors
                normx = norm(xk1-xk);
                normf = abs(fk1-fk);
                % display iterations of dca if actived
                if (obj.verbose == 1)
                    fprintf('%5d %19.5e %17.3e %13.3e\n',obj.iter,fk1,normx,normf);
                end
                % plotting if actived
                if (obj.plot==1)
                    myplotf(fk,fk1,obj.iter,obj.plotlinetype);
                end
                % check stopping
                if (normx < obj.tolx*(1+norm(xk1)) || normf < obj.tolf*(1+abs(fk1)))
                    if (obj.verbose == 1)
                        fprintf('------------------------------------------------------------\n');
                    end
                    obj.fopt = fk1;
                    obj.xopt = xk1;
                    status = setstatus(toc(cputime),obj.iter,0,'Successfully solved.');
                    return;
                end
                
                xk = xk1;
                fk = fk1;
                dfk = dfk1;
            end
            % maxiter exceed
            if (obj.verbose == 1)
                fprintf('------------------------------------------------------------\n');
            end
            obj.fopt = fk1;
            obj.xopt = xk1;
            status = setstatus(toc(cputime),obj.iter,1,'Max interation exceed.');
        end
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function myplotf(fk,fk1,iter,plotline)
    hold on;
    if iter == 1
        title('DCA iterations');
        xlabel('Iterations');
        ylabel('Objectives');
    end
    if iter>1
        plot([iter-1,iter], [fk,fk1],plotline);
        drawnow
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setting dca solution status
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status = setstatus(timer,iter,flag,info)
    status.time = timer;
    status.iter = iter;
    status.avgt = timer/iter;
    status.flag = flag;
    status.info = info;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% armijo
% classical armijo rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xx,ff] = armijo(obj,xk1,fk1,d,alpha,beta,sigma,eps)
    nd = norm(d);
    %ii=1;
    %idx= d<0;
    %alpha = min(-xk1(idx)./d(idx));
    while (alpha*nd>eps)
        xx=xk1+alpha*d;
        if checkfeas(xx)
            ff=obj.dcp.F.evalf(xx);
            delta=fk1 - ff - sigma*alpha^2*nd^2;
            if delta >0
                %fprintf('alpha:%.6f\n',alpha);
                return;
            end
        end
        alpha=beta*alpha;
        %ii = ii+1;
        %alpha = alpha / sqrt(ii);
    end
    xx=xk1;
    ff=fk1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adaptive armijo (method=0)
% increase many times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xx,ff,cur_alpha] = armijo_adaptive(obj,fk1,xk1,d,cur_alpha)
    % line search to get a better feasible solution
    nd = norm(d);
    alpha = max(min(cur_alpha,sqrt(2)/nd),1);
    beta = 0.3;
    while (alpha*nd>1e-5)
        xx=xk1+alpha*d;
        ff=obj.dcp.F.f(xx,1);
        delta=fk1 - ff; % - alpha^2*nd^2;
        if delta < 1e-8 || ~checkfeas(xx)
            alpha=beta*alpha; % back-tracking
        else
            cur_alpha=alpha;
            return;
        end
    end
    xx=xk1;
    ff=fk1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adaptive armijo (method=0)
% increase many times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xx,ff,cur_alpha] = armijo_adaptive_backforward(obj,xk1,d,cur_alpha)
    % line search to get a better feasible solution
    nd = norm(d);
    alpha = min(cur_alpha,sqrt(2)/nd);
    fk1=obj.dcp.F.f(xk1,1);
    beta = 0.3;
    backtrack=false;
    while (alpha*nd>1e-5)
        xx=xk1+alpha*d;
        ff=obj.dcp.F.f(xx,1);
        delta=fk1 - ff; % - alpha^2*nd^2;
        if delta < 1e-8 || ~checkfeas(xx)
            alpha=beta*alpha; % back-tracking
            backtrack=true;
        else
            if backtrack==true
                cur_alpha=alpha;
                return;
            else
                alpha=alpha/beta; % forward-tracking
            end
        end
    end
    xx=xk1;
    ff=fk1;
end

function [r,activeset]=checkfeas(x)
    %min(check(obj.dcp.C))>=-1e-7
    activeset = abs(x)<=1e-7;
    if abs(sum(x)-1) <=1e-7 && min(x) >= -1e-7
        r=true;
    else
        r=false;
    end
end