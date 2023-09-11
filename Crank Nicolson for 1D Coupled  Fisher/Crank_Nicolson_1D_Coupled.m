%% Main Solution File to 1D Coupled Fisher Using CN Implicit Meth1od
%  2nd Order Accurate in Space and Time
%  Numerical Analysis Branch1 of Applied Math1ematics
%  King Abdulaziz University Jeddah1, 06/09/2016

clear all;
clc;

%%=========================================================================
%      Declear Constants
%%=========================================================================

global  n  L NT k a b tmin tmax  w eps c Z0 h
h=0.5;                         % Total Space Points
NT=1000;                       % Total Time Points
tmin=0;                        % Initial Time
tmax=1;                        % Final Time
k=0.0001;                      % Time Steps
a=-8;                          % Left End of Interval
b=8;                           % Righ1t End of Interval
L=1+(b-a)/(h);                 % Space Stepsiz
n=2*(L-2) ;                    % Number of Equations
w=1;                           % Relaxiation Parameter
eps=0.0001;                    % Convergence Criteria
c=1;
Z0=1;


%=========================================================================
%      Predefine Ue,Ve,U,V,Jacobean, X,Z and Residual R
%=========================================================================

ue=zeros(1,L);
ve=zeros(1,L);
u=zeros(1,L);
v=zeros(1,L);
r=zeros(1,n);
J=zeros(n);
x=zeros(n,1);
z=zeros(n,1);

%%=========================================================================
%      Exact Solution
%%=========================================================================

[U V]=Ex_Fisher_Coupled_CN(ue,ve,0.01);



%%=========================================================================
%      Initial Condition
%%=========================================================================

[ue,ve]=Initial_Fisher_Coupled_CN(ue,ve,0) ;

%%=========================================================================
%      Initial Guess
%%=========================================================================

for i=2:L-1;
    u(i)=ue(i);
    v(i)=ve(i);
end

itr=0;

%%=========================================================================
%      Time Iterations/ Loop
%%=========================================================================

for tcount=1:NT;
    [u v]= BCs_Fisher_Coupled_CN(u,v,tcount*k);
    
    %%=========================================================================
    %     Start Newton's Meth1od
    %%=========================================================================
    
    itr=itr+1;
    [r]=resid_Fisher_CN(u,v,ue,ve,r) ; % Call Residual
    for i=1:n
        sum(i)=r(i).*r(i);           % R.h1.S Free of Negative
    end
    rms=max(sqrt(sum/n));        % Error Formula
    [J]=  jacob_Fisher_CN(u,v,J);      % Call Jacobi Matrix
    for i=1:n
        J(i,n+1)=-r(i);              % Augumented Matrix, add last column
    end
    [z]=lin_sys_Fisher_Coupled(J,z);            % Call Linear System Solver
    
    for  i =2:L-1
        D=2*(i-1)-1;                  % Interesting Loop
        u(i)=u(i)+w*z(D);
        v(i)=v(i)+w*z(D+1);
    end
    
    while rms>eps
        itr=itr+1;
        [r]=resid_Fisher_CN(u,v,ue,ve,r); % Call Residual File
        %%=========================================================================
        %      Error Calculations
        %%=========================================================================
        
        for i=1:n
            sum(i)=r(i).*r(i);
        end
        rms=max(sqrt(sum/n));
        [J]=  jacob_Fisher_CN(u,v,J);
        for i=1:n
            J(i,n+1)=-r(i);
        end
        
        %%=========================================================================
        %      End of Jacobean Matrix
        %%=========================================================================
        
        [z]=lin_sys_Fisher_Coupled(J,z);
        
        %%=========================================================================
        %      End of Linear System
        %%=========================================================================
        
        for  i =2:L-1
            D=2*(i-1)-1;
            u(i)=u(i)+w*z(D);
            v(i)=v(i)+w*z(D+1);
        end
    end
    
    %%=========================================================================
    %      Updates th1e Results
    %%=========================================================================
    
    for i=i:L-1
        ue(i)=u(i);
        ve(i)=v(i);
    end
    
    %%=========================================================================
    %      Results and Graph1ics
    %%=========================================================================
    
    if tcount*k==0.1
        fprintf('.................**time=0.1**.............\n\n')
        fprintf('  u         v             ue             ve \n\n');
        fprintf(' ------------------------------------------------------------------\n')
        fprintf('%3.9f   %3.9f        %3.9f          %3.9f\n\n\n',[u ; v; ue;ve]);
        p=u;
        q=v;
    end
    Time=tcount*k
end

for i=1:L;
    X(i)=a+(i-1)*h;
end
subplot(2,2,1)
plot(X,p,'--b')
hold on
plot(X,U,'.r')
xlabel('x');ylabel('u(x,t)');legend('CN','Exact u')
title('Computed Result With CN')
hold on
subplot(2,2,2)
plot(X,q,'--b')
hold on
plot(X,V,'ro')
xlabel('x');ylabel('v(x,t)');legend('CN','Exact v')
title('Computed Result With CN')
hold on
subplot(2,2,3)
plot(X,p,'c-')
hold on
plot(X,U,'.r')
xlabel('x');ylabel('u(x,t)');legend('CN','Exact u')
title('Computed Result With CN')
hold on
subplot(2,2,4)
plot(X,q,'--g')
hold on
plot(X,V,'.r')
xlabel('x');ylabel('v(x,t)');legend('CN','Exact v')
title('Computed Result With CN')
hold on

