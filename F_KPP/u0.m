function u = u0(x,cases)
if nargin == 1 
    cases = 1;
end
switch cases
    case 1
        v = 1/sqrt(2);
        u = (1 + exp(v*(x))).^(-1);
    case 2
        alpha = 6;
        u = (-1/2*tanh(alpha/(2*sqrt(2*alpha+4))*(x))+1/2).^(2/alpha);
    case 3
        a = 6;
        u = 1/2*(1+a)+(1/2-1/2*a)*tanh(sqrt(2)*(1-a)*x/4);
    otherwise
        u = zeros(length(x));
end
