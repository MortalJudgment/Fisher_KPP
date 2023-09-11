function dudt = odefn(t,u,N,dx,beta)

DD=zeros(N,1);
RR=zeros(N,1);

for i=1:N
    DD(i,1)=beta/(beta+u(i,1));
    RR(i,1)=beta*(1-u(i,1))*log(1 + u(i,1)/beta);
end

i=1;
dudt(i,1)=(DD(i+1,1).*(u(i+1,1)-u(i,1)))/dx^2+RR(i,1);
for i=2:N-1
    dudt(i,1)=((DD(i+1,1)+DD(i,1)).*(u(i+1,1)-u(i,1))-(DD(i-1,1)+DD(i,1)).*(u(i,1)-u(i-1,1)))/(2.*dx^2)+RR(i,1);
end
i=N;
dudt(i,1)=0.0;