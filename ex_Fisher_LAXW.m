%% Exact Solution to 1D Fisher Using Lax Weddroff Method
%  2nd Order Accurate in Space and Time
%  Produced by Shahid
%  Numerical Analysis Branch of Applied Mathematics
%  King Abdulaziz University Jeddah, 26/06/2016

function[ue]= ex_Fisher_LAXW(ue,t);

%==========================================================================
%      Declear Constants 
%==========================================================================

 global L a h alpha 
        
%==========================================================================
%      Declear and Finding of Exact Ue and Ve
%==========================================================================

for i=1:L;
    
    x(i)=a+(i-1)*h;
    ue(i)=1/(1+exp(sqrt(alpha/6)*x(i)-(5/6)*alpha*t))^2; 
end
return

%==========================================================================
%      Exact Solution Programme ended here.
%==========================================================================