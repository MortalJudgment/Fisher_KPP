clear all
close all
clc
beta = 1.0;
UU=1.0;
L=pi/(2*sqrt(1+beta));
dx=L/200;
N=round(L/dx+1);
u0=zeros(N,1);


for i=1:N
    x(i)=0.0+(i-1)*dx;
    u0(i) = beta*(exp(log(1 + UU/beta)*exp(-beta*0.0)*cos(sqrt(1+beta)*x(i)))-1.0);
    %u0(i)=1.0;
end
u0(N)=0.0;

tspan = [0 0.2 0.5 1];
[t,u] = ode45(@(t,u) odefn(t,u,N,dx,beta), tspan, u0);
figure


u_anal=zeros(length(tspan),N);

for tt=1:length(tspan)
    for i=1:N
        u_anal(tt,i)=beta*(exp(log(1 + UU/beta)*exp(-beta*tspan(tt))*cos(sqrt(1+beta)*x(i)))-1.0);
    end
end
plot(x,u0(:),'g--','LineWidth',2)
hold on
for i = 1:length(tspan)
    plot(x,u_anal(i,:),'k','LineWidth',2)
    plot(x,u(i,:),'g--','LineWidth',2)
end

xlabel('x')
ylabel('u(x,t)')
xlim([0 1.1*L])
ylim([0 UU])
xticks([0 L/2 L])
yticks([0 0.5 1.0 1.5 2.0])