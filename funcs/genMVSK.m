function MVSK=genMVSK(n,c,P,modtype)
    % generate MVSK model using different modeling language
    % MVSK=genMVSK(n,c,P)
    % n: nb of assets
    % c: investor's preferences
    % P: portfolio data
    % modtype: yalmip|polylab|sym|optimvar|sostools
    switch modtype
        case "yalmip"
            MVSK = proc_yalmip(n,c,P);
        case "polylab"
            MVSK = proc_polylab(n,c,P);
        case "sym"
            MVSK = proc_sym(n,c,P);
        case "optimvar"
            MVSK = proc_optimvar(n,c,P);
        case "sostools"
            MVSK = proc_sostools(n,c,P);
        otherwise
            MVSK = proc_sym(n,c,P);
    end
end

function MVSK = proc_yalmip(n,c,P)
    MVSK.x=sdpvar(n,1); % decision variables
    MVSK.Cons=[sum(MVSK.x)==1;MVSK.x>=0]; % simplex constraint
    % generate datas
    MVSK.Data.P=P; % portfolio
    MVSK.Data.IDXS=genskewnessidx_s(P); % index sets for S
    MVSK.Data.IDXK=genkurtosisidx_s(P); % index sets for K
    MVSK.Data.c=c; % investor's preferences
    % generate functions
    MVSK.m1=m1(MVSK.x,P.mu);
    MVSK.m2=m2(MVSK.x,P.Sigma);
    MVSK.m3=m3_ss(MVSK.x,P);
    MVSK.m4=m4_ss(MVSK.x,P);
    MVSK.fobj=fobj(MVSK.x,c,MVSK.m1,MVSK.m2,MVSK.m3,MVSK.m4);
end

function MVSK = proc_optimvar(n,c,P)
    MVSK.x=optimvar('x',n,1,'LowerBound',0,'UpperBound',1);
    MVSK.Cons= sum(MVSK.x)==1; % simplex constraint
    % generate datas
    MVSK.Data.P=P; % portfolio
    MVSK.Data.IDXS=genskewnessidx_s(P); % index sets for S
    MVSK.Data.IDXK=genkurtosisidx_s(P); % index sets for K
    MVSK.Data.c=c; % investor's preferences
    % generate functions
    MVSK.m1=m1(MVSK.x,P.mu);
    MVSK.m2=m2(MVSK.x,P.Sigma);
%     MVSK.m3=m3_ss(MVSK.x,P);
    MVSK.m4=m4_ss(MVSK.x,P);
    MVSK.fobj=fobj(MVSK.x,c,MVSK.m1,MVSK.m2,MVSK.m3,MVSK.m4);
end

function MVSK = proc_sym(n,c,P)
    MVSK.x=sym('x',[n,1]); % decision variables
    MVSK.Cons=[sum(MVSK.x)==1;MVSK.x>=0]; % simplex constraint
    % generate datas
    MVSK.Data.P=P; % portfolio
    MVSK.Data.IDXS=genskewnessidx_s(P); % index sets for S
    MVSK.Data.IDXK=genkurtosisidx_s(P); % index sets for K
    MVSK.Data.c=c; % investor's preferences
    % generate functions
    MVSK.m1=expand(m1(MVSK.x,P.mu));
    MVSK.m2=expand(m2(MVSK.x,P.Sigma));
    MVSK.m3=expand(m3_ss(MVSK.x,P));
    MVSK.m4=expand(m4_ss(MVSK.x,P));
    MVSK.fobj=expand(fobj(MVSK.x,c,MVSK.m1,MVSK.m2,MVSK.m3,MVSK.m4));
end

function MVSK = proc_sostools(n,c,P)
    MVSK.x=mpvar('x',[n,1]); % decision variables
    MVSK.Cons.Aeq = ones(1,n);
    MVSK.Cons.beq = 1;
    MVSK.Cons.lb = zeros(n,1);
    MVSK.Cons.ub = ones(n,1);
    % generate datas
    MVSK.Data.P=P; % portfolio
    MVSK.Data.IDXS=genskewnessidx_s(P); % index sets for S
    MVSK.Data.IDXK=genkurtosisidx_s(P); % index sets for K
    MVSK.Data.c=c; % investor's preferences
    % generate functions
    MVSK.m1=m1(MVSK.x,P.mu);
    MVSK.m2=m2(MVSK.x,P.Sigma);
    MVSK.m3=m3_ss(MVSK.x,P);
    MVSK.m4=m4_ss(MVSK.x,P);
    MVSK.fobj=fobj(MVSK.x,c,MVSK.m1,MVSK.m2,MVSK.m3,MVSK.m4);
end

function MVSK = proc_polylab(n,c,P)
    MVSK.x=MPOLY.mpolyvars(n); % decision variables
    % define constraints
    MVSK.Cons.Aeq = ones(1,n);
    MVSK.Cons.beq = 1;
    MVSK.Cons.lb = zeros(n,1);
    MVSK.Cons.ub = ones(n,1);
    % generate datas
    MVSK.Data.P=P; % portfolio
    MVSK.Data.IDXS=genskewnessidx_s(P); % index sets for S
    MVSK.Data.IDXK=genkurtosisidx_s(P); % index sets for K
    MVSK.Data.c=c; % investor's preferences
    % generate functions
    if c(1)~=0
        MVSK.m1=m1(MVSK.x,P.mu);
    else
        MVSK.m1=MPOLY.zeros(n,1,1);
    end
    if c(2)~=0
        MVSK.m2=m2(MVSK.x,P.Sigma);
    else
        MVSK.m2=MPOLY.zeros(n,1,1);
    end
    if c(3)~=0
        MVSK.m3=m3_ss(MVSK.x,P);
    else
        MVSK.m3=MPOLY.zeros(n,1,1);
    end
    if c(4)~=0
        MVSK.m4=m4_ss(MVSK.x,P);
    else
        MVSK.m4=MPOLY.zeros(n,1,1);
    end
    MVSK.fobj=fobj(MVSK.x,c,MVSK.m1,MVSK.m2,MVSK.m3,MVSK.m4);
end

