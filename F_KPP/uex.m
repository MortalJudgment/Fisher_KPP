function u = uex(x,t,cases)
if nargin == 2
    cases = 1;
end

switch cases
    case 1
        v = 1/sqrt(2);
        u = (1 + exp(v*(x-v*t))).^(-1);
    case 2
        alpha = 1;
        u = (-1/2*tanh(alpha/(2*sqrt(2*alpha+4))*(x-(alpha+4)/(sqrt(2*alpha+4))*t))+1/2).^(2/alpha);
    case 3
        a = 1/2;
        u = 1/2*(1+a)+(1/2-1/2*a)*tanh(sqrt(2)*(1-a)*x/4 + (1-a^2)/4*t);
    otherwise
        u = zeros(length(x));
end
