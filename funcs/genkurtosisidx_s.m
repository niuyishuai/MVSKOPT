function IDX = genkurtosisidx_s(P)
% generate index sets for K
n=P.n;
tolzero=1e-8;
N=1:n;
Ps=[];
Ph=[];
Qh=[];
R=[];
IPK=[];
IMK=[];
JPK=[];
JMK=[];
JHPK=[];
JHMK=[];
KHPK=[];
KHMK=[];
LPK=[];
LMK=[];
for i=N % in N
    if (cokurtosis(i,i,i,i,P)>tolzero)
        IPK=[IPK i];
    elseif (cokurtosis(i,i,i,i,P)<-tolzero)
        IMK=[IMK i];
    end
    for j=N
        if(j~=i) % in P
            Ps=[Ps [i;j]];
            if (cokurtosis(i,i,i,j,P)>tolzero)
                JPK=[JPK [i;j]];
            elseif (cokurtosis(i,i,i,j,P)<-tolzero)
                JMK=[JMK [i;j]];
            end
        end
        if(j>i) % in Ph
            Ph=[Ph [i;j]];
            if (cokurtosis(i,i,j,j,P)>tolzero)
                JHPK=[JHPK [i;j]];
            elseif (cokurtosis(i,i,j,j,P)<-tolzero)
                JHMK=[JHMK [i;j]];
            end
        end
        for k=N
            if (j<k && i~=j && i~=k) % in Qh
                Qh=[Qh [i;j;k]];
                if (cokurtosis(i,i,j,k,P)>tolzero)
                    KHPK=[KHPK [i;j;k]];
                elseif (cokurtosis(i,i,j,k,P)<-tolzero)
                    KHMK=[KHMK [i;j;k]];
                end
            end
            for l=N
                if (i<j && j<k && k<l) % in R
                    R=[R [i;j;k;l]];
                    if (cokurtosis(i,j,k,l,P)>tolzero)
                        LPK=[LPK [i;j;k;l]];
                    elseif (cokurtosis(i,j,k,l,P)<-tolzero)
                        LMK=[LMK [i;j;k;l]];
                    end
                end
            end
        end
    end
end
IDX.N=N;
IDX.P=Ps;
IDX.Ph=Ph;
IDX.Qh=Qh;
IDX.R=R;
IDX.IPK=IPK;
IDX.IMK=IMK;
IDX.JPK=JPK;
IDX.JMK=JMK;
IDX.JHPK=JHPK;
IDX.JHMK=JHMK;
IDX.KHPK=KHPK;
IDX.KHMK=KHMK;
IDX.LPK=LPK;
IDX.LMK=LMK;
end