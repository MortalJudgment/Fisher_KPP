% Finite Difference method for Fisher-KPP equation on 1D
% Using explicit Forward in Time Central in Space method (FTCS)
% Consider on Domain Omega = (0 1), t>0
% with initial condition u0(x)
% and Dirichlet boudary condition u(0,t)=g1(t); u(1,t)=g2(t)
%-------------------------------------------------------------------------%
%----------------------- Fisher-KPP Equation: ----------------------------%
%------         u_t = D*u_xx + k*u(1-u)     ,in [a,b]              -------%
%-------------------------------------------------------------------------%
clc
clear all
close all
%----------------%
% Domain of x
ax = -10.0;
bx = 20.0;
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
u1 = (u0(x));
u_dis = zeros(1,length(x));
% theta = 1/2;
%------------------------------------------------------------------------%
figh = figure;
for jj=1:200
    clf
    t_k = jj*dt;
%     %---------- Dicrete solution -----------%
%     
    for ii = 2:Nx
        u_dis(ii) = u1(ii) + dt/(dx^2)*(u1(ii+1) - 2*u1(ii) + u1(ii-1)) + dt*u1(ii)*(1 - u1(ii));
    end
    u_dis(1) = uex(ax,t_k);
    u_dis(Nx+1) = uex(bx,t_k);
    u1 = u_dis;
    %----------- Exact solution ------------%
    u_ex = zeros(Nx+1,1);
    for ii=1:Nx+1
        u_ex(ii) = uex(x(ii),t_k);
    end
    %---------------------------------------%
    %-------------- Drawing ----------------%
    %---------------------------------------%
    hold on
    plot(x,u_dis,'o r','LineWidth',2,'MarkerSize',4);
    plot(x,u_ex,'o b');
    
%     legend('discrete solution','exact solution')
%     axis([ax bx 0 1])
    title(['Particle at t = ', num2str(jj),'*dt seconds'])
    pause(10^(-2))
    movieVector(jj) = getframe(figh, [10 10 520 400]);
end

myWriter = VideoWriter('Fisher-KPP Equation Compare','MPEG-4');
myWriter.FrameRate = 20;
open(myWriter);
writeVideo(myWriter,movieVector);
close(myWriter);
hold off