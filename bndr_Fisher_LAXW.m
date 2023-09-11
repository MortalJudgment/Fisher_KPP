%% Boundary Conditions to 1D Fisher Using Lax Wendroff Method
%  2nd Order Accurate in Space and Time
%  Produced by Shahid
%  Numerical Analysis Branch of Applied Mathematics
%  King Abdulaziz University Jeddah, 26/06/2016

function[ue]= bndr_Fisher_LAXW(ue,t);

%%=========================================================================
%      Declear Constants 
%%=========================================================================

global L a b alpha

%%=========================================================================
%      Declear Boundary values
%%=========================================================================

x(1)=a;

       ue(1)=1/(1+exp(sqrt(alpha/6)*x(1)-(5/6)*alpha*t))^2;
       
x(L)=b;

      ue(L)=1/(1+exp(sqrt(alpha/6)*x(L)-(5/6)*alpha*t))^2;

return

%%=========================================================================
%      Boundary Conditions Programme ended here.
%%=========================================================================
