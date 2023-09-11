% Finite Difference method for Fisher-KPP equation on 1D
% Using explicit Forward in Time Central in Space method (FTCS)
% Consider on Domain Omega = (0 1), t>0
% with initial condition u0(x)
% and Dirichlet boudary condition u(0,t)=g1(t); u(1,t)=g2(t)
%-------------------------------------------------------------------------%
%----------------------- Fisher-KPP Equation: ----------------------------%
%--------         u_t = u_xx + u(1-u)     ,in [a,b]              ---------%
%-------------------------------------------------------------------------%
clc
clear all
close all
%----------------%
% Domain of x
ax = -10.0;
bx = 10.0;
%Coefficience
beta = 1/16;
%----------------%
Nx = 100 ;
%-------------------------------------------------------------------------%

dx = (bx-ax)/Nx;
dt = 0.01;
% Create the mesh point
x = zeros(Nx+1,1);
for ii=1:Nx+1
    x(ii) = ax+(ii-1)*dx; 
end
%-----------------%
%Initial condition
u1 = (u0(x,2));
u_dis = zeros(length(x),1);
u_star = zeros(length(x),1);
% theta = 1/2;
%------------------------------------------------------------------------%
for jj=1:300
    clf
    t_k = jj*dt;
    %---------- Dicrete solution -----------%
    alpha=1;
    for ii = 2:Nx
        u_star(ii) = 0.5*(u1(ii+1) + u1(ii-1)) - 0.5*dt/dx^2*(u1(ii+1) - 2*u1(ii) + u1(ii-1)) ...
            + 0.5*alpha*dt*u1(ii)*(1-u1(ii));
    end
    u_star(1) = uex(ax,t_k,2);
    u_star(Nx+1) = uex(bx,t_k,2);
    for kk = 2:Nx
        u_dis(kk) = u_star(kk) + dt/dx^2*(u_star(kk+1) - 2*u_star(kk) + u_star(kk-1)) ...
            + 0.5*alpha*dt*u_star(kk)*(1 - u_star(kk));
    end
    u_dis(1) = uex(ax,t_k,2);
    u_dis(Nx+1) = uex(bx,t_k,2);
    
    u1 = u_dis;
    %----------- Exact solution ------------%
    u_ex = zeros(Nx+1,1);
    for ii=1:Nx+1
        u_ex(ii) = uex(x(ii),t_k,2);
    end
    %---------------------------------------%
    %-------------- Drawing ----------------%
    %---------------------------------------%
    hold on
    plot(x,u_dis,'o r','LineWidth',2,'MarkerSize',4);
    plot(x,u_ex,'b');
    
    %     legend('discrete solution','exact solution')
    % axis([0 1 -1 1])
    title(['Particle at t = ', num2str(jj),'*dt seconds'])
    pause(10^(-2))
end
hold off