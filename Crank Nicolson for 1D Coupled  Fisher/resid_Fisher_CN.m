%% Residuals of 1D Coupled Fisher Using CN Implicit Method
%  2nd Order Accurate in Space and Time
%  Numerical Analysis Branch of Applied Mathematics
%  King Abdulaziz University Jeddah, 06/09/2016

function[r]= resid_Fisher_CN(u,v,u1,v1,r);

%==========================================================================
%      Declear Constants
%==========================================================================

global L k h

%==========================================================================
%      Declear Residual Equations.........AX=B    ==>   r=B-AX
%==========================================================================

for i=2:L-1;
    
    r((2*i)-3)=((u(i)-u1(i))/k)-0.5*((u(i-1)-2*u(i)+u(i+1))/(h*h)+(u1(i-1)-2*u1(i)+u1(i+1))/(h*h))-...
        0.5*((u(i)+ u1(i)))*(1-0.5*((u(i)+u1(i))))*(1-0.5*((v(i)+v1(i))));
    
    r((2*i)-2)=((v(i)-v1(i))/k)-0.5*((v(i-1)-2*v(i)+v(i+1))/(h*h)+(v1(i-1)-2*v1(i)+v1(i+1))/(h*h))-...
        0.5*((v(i)+ v1(i)))*(1-0.5*((v(i)+v1(i))))*(1-0.5*((u(i)+u1(i))));
    
end
end

%==========================================================================
%      Residual Programme ended here.
%==========================================================================


