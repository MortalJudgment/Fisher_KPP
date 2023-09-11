%% Linear System of 1D Coupled Fisher Using CN Implicit Method
%  2nd Order Accurate in Space and Time
%  Numerical Analysis Branch of Applied Mathematics
%  King Abdulaziz University Jeddah, 07/09/2016

function[y]= lin_sys_Fisher_Coupled(A,y)

%==========================================================================
%      Declear Constants
%==========================================================================

global    n

%==========================================================================
%      Declear Linear System and Its Solution Solver.........AX=B
%==========================================================================

for i=1:n-1
    t=abs(A(i,i));
    p=i;
    for j=i+1:n
        if abs(A(j,i))>t
            p=j;
            t=abs(A(j,i));
        end
    end
    if p~=i
        for j=i:n+1
            temp=A(i,j);
            A(i,j)=A(p,j);
            A(p,j)=temp;
        end
    end
    for j=i+1:n
        t1=(A(j,i))/(A(i,i));
        if abs(t1)<=0
            t1=0;
        else
            for k=i:n+1
                A(j,k)=A(j,k)-t1.*A(i,k);
            end
        end
        
    end
end
if abs(A(n,n))<=0
    fprintf('Unique Solution Does Not Exist');
else
    %//using back substitution method
    
    y(n)=A(n,n+1)/A(n,n);
    for i=1:n-1
        j=n-i;
        for k=j+1:n
            A(j,n+1)=A(j,n+1)-A(j,k)*y(k);
            y(j)=(A(j,n+1))/(A(j,j));
        end
    end
end

%==========================================================================
%      Linear System Solver Programme ended here.
%==========================================================================
