function f=m4_ss(x,P)
    % create m4 (portfolio kurtosis)
    % f=m4_ss(x,P)
    n=P.n;
    v1=0;
    v2=0;
    v2b=0;
    v3=0;
    v4=0;
    %figure(1)
    %hold on
    %ss=0;
    %tic
    for i=1:n
        v1=v1+cokurtosis(i,i,i,i,P)*x(i)^4;
        %   mydraw;
        for j=1:n
            if(j~=i)
                v2=v2+cokurtosis(i,i,i,j,P)*x(i)^3*x(j);
                %          mydraw;
            end
            if(j>i)
                v2b=v2b+cokurtosis(i,i,j,j,P)*x(i)^2*x(j)^2;
                %         mydraw;
                for k=j+1:n
                        v3=v3+cokurtosis(i,i,j,k,P)*x(i)^2*x(j)*x(k);
                        v3=v3+cokurtosis(i,j,j,k,P)*x(i)*x(j)^2*x(k);
                        v3=v3+cokurtosis(i,j,k,k,P)*x(i)*x(j)*x(k)^2;                        
                        %            mydraw;
                    for l=k+1:n
                            v4=v4+cokurtosis(i,j,k,l,P)*x(i)*x(j)*x(k)*x(l);
                            %               mydraw;
                    end
                end
                
            end
        end
    end
    f=v1+4*v2+6*v2b+12*v3+24*v4;
    %mydraw;
    %    function mydraw
    %       ss=ss+1;
    %       plot(ss,toc,'r-o');
    %       drawnow;
    %    end
end