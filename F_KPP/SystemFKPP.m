% Finite Difference method for Fisher-KPP equation on 1D
% Using explicit Forward in Time Central in Space method (FTCS)
% Consider on Domain Omega = (0 1), t>0
% with initial condition u0(x)
% and Dirichlet boudary condition u(0,t)=g1(t); u(1,t)=g2(t)
%-------------------------------------------------------------------------%
%----------------------- Fisher-KPP Equation: ----------------------------%
%---                a_t = a_xx + r*a(1-a-i)                            ---%
%---                i_t = a + r*a(a+i)                                 ---%                               
%-------------------------------------------------------------------------%
clc
clear all
close all
%----------------%
% Domain of x
ax = 0.0;
bx = 1000.0;
Nx = 1000 ;
% Step size
dx = (bx-ax)/Nx;
dt = 0.0005;
r = 1;
% Create the mesh point
x = zeros(Nx+1,1);
for ii=1:Nx+1
    x(ii) = ax + (ii-1)*dx; 
end
%-----------------%
%Initial condition 
u1 = zeros(1,length(x));
v1 = zeros(1,length(x));
for ii=1:Nx+1
    u1(ii) = u0(x(ii)); 
end
for ii=1:Nx+1
    v1(ii) = v0(x(ii)); 
end
u_dis = zeros(1,length(x));
v_dis = zeros(1,length(x));
%------------------------------------------------------------------------%
figh = figure;
for jj=1:10000
    clf
    t_k = jj*dt;
    %---------- Dicrete solution -----------%
    
    for ii = 2:Nx
        u_dis(ii) = u1(ii) + dt/dx^2*(u1(ii+1) - 2*u1(ii)+u1(ii-1)) + (r*dt)*u1(ii)*(1 - u1(ii) - v1(ii));
    end
    u_dis(1) = u_dis(2);
    u_dis(Nx+1) = 0;
    for ii = 2:Nx
        v_dis(ii) = v1(ii) + dt*u1(ii) + dt*r*u1(ii)*(u1(ii) + v1(ii));
    end
    v_dis(1) = v_dis(2);
    v_dis(Nx+1) = v_dis(Nx);
    
%     u1 = u_dis;
%     v1 = v_dis;
    
    u1(2:Nx+1) = u_dis(1:Nx);
    u1(1) = u1(2);
    v1(2:Nx+1) = v_dis(1:Nx);
    v1(1) = v1(2);
    
    for ii=1:Nx+1
        u_dis(ii) = min(u_dis(ii));
        v_dis(ii) = min(1 - u_dis(ii),v_dis(ii));
    end
    
%     %----------- Exact solution ------------%
%     u_ex = zeros(Nx+1,1);
%     for ii=1:Nx+1
%         u_ex(ii) = uex(x(ii),t_k);
%     end
%     v_ex = zeros(Nx+1,1);
%     for ii=1:Nx+1
%         v_ex(ii) = vex(x(ii),t_k);
%     end
    %---------------------------------------%
    %-------------- Drawing ----------------%
    %---------------------------------------%
    hold on
    plot(x(2:Nx),u_dis(2:Nx),'b','LineWidth',2,'MarkerSize',4);
    plot(x(2:Nx),v_dis(2:Nx),'r','LineWidth',2,'MarkerSize',4);
%     plot(x,u_ex,'b');
%     plot(x,v_ex,'g');
    
%     pause(0.005);
%     legend('discrete solution','exact solution')
%     axis([ax bx 0 1])
    title(['Particle at t = ', num2str(jj),'*dt seconds'])
    pause(10^(-2))
    movieVector(jj) = getframe(figh, [10 10 520 400]);
end

myWriter = VideoWriter('system of Fisher-KPP Equation','MPEG-4');
myWriter.FrameRate = 20;
open(myWriter);
writeVideo(myWriter,movieVector);
close(myWriter);
hold off