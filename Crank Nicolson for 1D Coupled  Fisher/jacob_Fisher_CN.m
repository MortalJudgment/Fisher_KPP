%% Jacobi Matrix of 1D Coupled Fisher Using CN Implicit Method
%  Produced by Shahid
%  Numerical Analysis Branch of Applied Mathematics
%  King Abdulaziz University Jeddah, 06/09/2016

function[J]=jacob_Fisher_CN( u, v ,J)

%==========================================================================
%      Declear Constants
%==========================================================================

global L k  h n

%==========================================================================
%      Declear and Finding of Jacobi Matrix
%==========================================================================

for i=2:L-1
    
    D= 2*(i-1)-1;
    
    J(D,D)= (1/k) + 1/(h*h)- 0.5*(1-0.5*(u(i)))*(1-0.5*(v(i)))+0.25*(1-0.5*(v(i)))*(u(i)); %//deri w.r.t. u[i] of 1st residual
    J(D,D+1)    = 0.25*(1-0.5*(u(i)))*(u(i));   %//deri w.r.t. v[i]
    
    D1 = D+2;
    if D1<=n
        J(D,D1)     = -1/(2*h*h);  %//deri w.r.t. u[i+1]
    end
    
    D1 = D-2;
    if D1>0
        J(D,D1)     = 0;   %//deri w.r.t. v[i+1]
    end
    D1 = D+2*L-4;
    if(D1<=n)
        J(D,D1)     = -1/(2*h*h);  %//deri w.r.t. u[i-1]
    end
    D1 = D-2*L+4;
    if(D1>0)
        J(D,D1)     =0;     %//deri w.r.t. v[i-1]
    end
    
    %====================2nd Block of Jacobi=================================
    
    
    J(D+1,D)    =  0.25*(1-0.5*(v(i)))*(v(i));  %//deri w.r.t. u[i] of 2nd residual
    J(D+1,D+1)  = (1/k) + 1/(h*h)- 0.5*(1-0.5*(u(i)))*(1-0.5*(v(i)))+0.25*(1-0.5*(u(i)))*(v(i));   %//deri w.r.t. v[i]
    
    D1 = D+2;
    if(D1<=n)
        J(D+1,D1+1) = 0;    %//deri w.r.t. u[i+1]
    end
    D1 = D-2;
    if(D1>0)
        J(D+1,D1+1) = -1/(2*h*h);%//deri w.r.t. v[i+1]
    end
    D1 = D+2*L-4;
    if(D1<=n)
        J(D+1,D1+1) =0;     %//deri w.r.t. u[i-1]
    end
    D1 = D-2*L+4;
    if(D1>0)
        J(D+1,D1+1) =-1/(2*h*h);    %//deri w.r.t. v[i-1]
    end
end

%==========================================================================
%      Jacobi Programme ended here.
%==========================================================================
