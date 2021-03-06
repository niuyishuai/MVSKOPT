%clear;
n_min=100;% min number of assets
n_max=100;% max number of assets
T=30;% number of periods
nbprobsforeachn=1;% number of problems for each n
range_n=n_min:2:n_max; % range of n
range_i=1:nbprobsforeachn; % range of i for each n
%% Generate MVSK models' data (randomly)
if ~exist('datas','dir')==0
   mkdir('datas');
end
for n=range_n
    for i=range_i
        tic
        fprintf('Creating co-moments for (n,T,i) = (%d,%d,%d) ... \n',n,T,i);
        P=genRandPortfolio_mo(n,T); % generate random portfolio data
        switch i % create investor's prefs
            case 1
                c=[10,1,10,1];
            case 2
                c=[1,10,1,10];
            case 3
                c=[10,10,10,10];
        end
        x0=ones(n,1)/n; % use 1/n for initialization.
        fname=sprintf('%d_%d_%d',n,T,i);
        fprintf('Saving portfolio data to file ...\n');
        save(['datas//',fname,'.mat'],'n','T','c','P','x0');
        toc
    end
end
fprintf("End of all portfolio data and models generation.\n");

%% initialize some parameters for tests
% flag setting
testudca=1;
testubdca=1;

% tolerance
tolf=1e-6;
tolx=sqrt(tolf);

% initialize output lists
nbprobs=numel(range_n)*numel(range_i);
zerolst=zeros(1,nbprobs);
genprob_times=zerolst;
if testudca
    udca_times=zerolst;
    udca_objs=zerolst;
    udca_iters=zerolst;
end
if testubdca
    ubdca_times=zerolst;
    ubdca_objs=zerolst;
    ubdca_iters=zerolst;
end

%% solve all MVSK models using UDCA and UBDCA
counter = 1;
for n=range_n
    for i=range_i
        % read data
        fname=sprintf('%d_%d_%d',n,T,i);
        %fprintf('Loading MVSK model from file %s.mat\n',fname);
        %load(['datas//',fname,'.mat']); % load n, T, c, P, x0
        %fprintf('End MVSK model loading.\n');
        
        % Initialize DCAM object for portfolio model
        % compute rho
        fprintf('Initializing DC programming model.\n');
        rho = 2*c(2)*norm(P.Sigma,'inf')+6*c(3)*max(norm(P.coskewness,'inf')) ...
            + 12*c(4)*max(norm(P.cokurtosis,'inf'));
        
        % create a dc function object
        dcf=dcfuncpoly;
        dcf.f=@(x,opt)fobj_eval(x,c,P,opt);
        
        % create a dc problem object
        mydcp=dcppoly(dcf,P.Cons);
        fprintf('End of initialization.\n');
        
        %%
        % Test UDCA solver
        if testudca
            fprintf('Solving MVSK model using UDCA.\n');
            % create and initialize a dca object
            mydca = dcapoly(mydcp,x0);
            mydca.rho=rho;
            mydca.P=P;
            mydca.c=c;
            mydca.plot=0;
            mydca.tolf=tolf;
            mydca.tolx=tolx;
            mydca.verbose = 0;
            mydca.convexsolver='bpppa';
            mydca.linesearch=0;
            
            % solve model using dca
            status=mydca.optimize();
  
            % get results
            udcatime=status.time;
            udca_times(counter) = udcatime;
            udca_objs(counter) = mydca.fopt;
            udca_iters(counter) = status.iter;    
            fprintf('Solution for UDCA: time %.3f sec, obj %.4e iters %d\n',udca_times(counter),udca_objs(counter),status.iter);
        end
        
        %%
        % Test UBDCA solver
        if testubdca
            fprintf('Solving MVSK model using UBDCA.\n');
            % create a dca object
            mydca = dcapoly(mydcp,x0);
            mydca.rho=rho;
            mydca.P=P;
            mydca.c=c;
            mydca.plot=0;
            mydca.tolf=tolf;
            mydca.tolx=tolx;
            mydca.verbose = 0;
            mydca.convexsolver='bpppa';
            mydca.linesearch=1;
            
            % solve model using ubdca
            status=mydca.optimize();
            
            % get results
            ubdcatime=status.time;
            ubdca_times(counter) = ubdcatime;
            ubdca_objs(counter) = mydca.fopt;
            ubdca_iters(counter) = status.iter;
            
            fprintf('Solution for UBDCA: time %.3f sec, obj %.4e iters %d\n',ubdca_times(counter),ubdca_objs(counter),status.iter);
        end
        counter = counter + 1;
    end
end

fprintf('All test finished!\n');