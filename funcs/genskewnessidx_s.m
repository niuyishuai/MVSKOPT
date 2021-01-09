function IDX = genskewnessidx_s(P)
% generate index sets for S
n=P.n;
tolzero=1e-8;
N=1:n;
Ps=[];
Q=[];
IPS=[];
IMS=[];
JPS=[];
JMS=[];
KPS=[];
KMS=[];
for i=N % in set N
    if (coskewness(i,i,i,P)>tolzero)
        IPS=[IPS i];
    elseif (coskewness(i,i,i,P)<-tolzero)
        IMS=[IMS i];
    end
    for j=N
        if(j~=i) % in set P
            Ps=[Ps [i;j]];
            if (coskewness(i,i,j,P)>tolzero)
                JPS=[JPS [i;j]];
            elseif (coskewness(i,i,j,P)<-tolzero)
                JMS=[JMS [i;j]];
            end
        end
        for k=N
            if (i<j && j<k) % in set Q
                Q=[Q [i;j;k]];
                if (coskewness(i,j,k,P)>tolzero)
                    KPS=[KPS [i;j;k]];
                elseif (coskewness(i,j,k,P)<-tolzero)
                    KMS=[KMS [i;j;k]];
                end
            end
        end
    end
end
IDX.N=N;
IDX.P=Ps;
IDX.Q=Q;
IDX.IPS=IPS;
IDX.IMS=IMS;
IDX.JPS=JPS;
IDX.JMS=JMS;
IDX.KPS=KPS;
IDX.KMS=KMS;
end