%% Main Solution File to 1D Fisher Using Lax Wendroff Method
%  2nd Order Accurate in Space and Time
%  Produced by Shahid
%  Numerical Analysis Branch of Applied Mathematics
%  King Abdulaziz University Jeddah, 26/06/2016

clear all;
clc;

%%=========================================================================
%      Declear Constants
%%=========================================================================

global L NT k a b tmin tmax h alpha R_1

L=5;                          % Total Space Points
NT=1000;                    % Total Time Points
tmin=0;                        % Initial Time
tmax=1;                        % Final Time
k=0.0001         % Time Steps
a=-8;                          % Left End of Interval
b=8;                           % Right End of Interval
h=(b-a)/(L-1);                 % Space Stepsize
eps=0.0001;                 % Convergence Criteria
alpha=1;                     % Reaction Constant
R_1=k/(h*h);

%%=========================================================================
%      Stability Condition
%%=========================================================================

if R_1 > 1
    fprintf('Stability does not Satisfied \n');
else
    
    %%=========================================================================
    %      Initialize Variables in Codeing
    %%=========================================================================
    
    ue=zeros(1,L);
    u=zeros(1,L);
    Ue=zeros(1,L);
    
    %%=========================================================================
    %      Exact and Initial
    %%=========================================================================
    
    [U]=ex_Fisher_LAXW(ue,0.001);
    
    [ue]=int_Fisher_LAXW(ue,0.0);
    
    for i=1:L
        X(i)=a+(i-1)*h;
    end
    itr=0;
    %%=========================================================================
    %     Time Loop
    %%=========================================================================
    
    for tcount=1:NT
        [u]= bndr_Fisher_LAXW(u,tcount*k);
        
        %%=========================================================================
        %      Discritization
        %%=========================================================================
        
        for i=2:L-1;
            Ue(i)=0.5*(ue(i+1)+ue(i-1))+0.5*R_1*(ue(i+1)-2*ue(i)+ue(i-1))+0.5*k*alpha*ue(i)*(1-ue(i));
            u(i)=ue(i)+R_1*(Ue(i+1)-2*Ue(i)+Ue(i-1))+alpha*k*Ue(i)*(1-Ue(i));
        end
        
        %%=========================================================================
        %      Error Calculations
        %%=========================================================================
        
        for i= 1:L
            Error1=abs(u(i)-ue(i));
            Error=max(Error1);
        end
        
        %%=========================================================================
        %      Convergance Criteria
        %%=========================================================================
        
        if Error <= eps;
            fprintf('Successful at Iteration = %g \n ',itr);
            break
        else
            
            %%=========================================================================
            %      Update Output After Convergance
            %%=========================================================================
            
            for i=1:L-1
                ue(i)=u(i);
            end
            itr=itr+1;
        end
    end
    Time=tcount*k
end

%%=========================================================================
%      Results and Graphics
%%=========================================================================
subplot(2,2,1);
plot(X,u,'b.')
hold on
plot(X,U,'r.')
xlabel('x');
ylabel('u(x,t)');
legend('u_{app.}', 'u_{Exact}');
title('Computed Solution With Crank Nicolson')
hold on
subplot(2,2,2);
plot(X,u,'b.');
hold on
plot(X,U,'r.');
xlabel('x');
ylabel('u(x,t)');
legend('u_{app.}', 'u_{Exact}');
title('Computed Solution With Crank Nicolson')
hold on
xlabel('x');ylabel('u(x,t)');legend('Lax Wendroff','Exact u','Lax Wendroff','Exact u','Lax Wendroff','Exact u')
title('Computed Result With Lax Wendroff')
hold on
subplot(2,2,3);
plot(X,u,'g');
hold on
plot(X,U,'r');
xlabel('x');
ylabel('u(x,t)');
legend('u_{app.}', 'u_{Exact}');
title('Computed Solution With Crank Nicolson')
hold on
subplot(2,2,4);
plot(X,u,'g');
hold on
plot(X,U,'r');
xlabel('x');
ylabel('u(x,t)');
legend('u_{app.}', 'u_{Exact}');
title('Computed Solution With Crank Nicolson')
hold on


%%=========================================================================
%      Lax Wendroff Programme ended here
%%=========================================================================





