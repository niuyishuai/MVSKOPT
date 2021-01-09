function f=m3_ss(x,P)
    % create m3 (portfolio skewness)
    % f=m3_ss(x,P)
    n=P.n;
    v1=0;
    v2=0;
    v3=0;
    %figure(1);
    %ss=0;
    %hold on;
    %tic
    for i=1:n
        v1=v1+coskewness(i,i,i,P)*x(i)^3;
        %    mydraw;
        for j=1:n
            if(j~=i)
                v2=v2+coskewness(i,i,j,P)*x(i)^2*x(j);
                %            mydraw;
            end
            if (i<j)
                for k=j+1:n
                    v3=v3+coskewness(i,j,k,P)*x(i)*x(j)*x(k);
                    %                mydraw;
                end
            end
        end
    end
    f=v1+3*v2+6*v3;
    %mydraw;
    %    function mydraw
    %       ss=ss+1;
    %       plot(ss,toc,'r-o');
    %       drawnow;
    %    end
end
