function [P, T, rho] = ISA(h)
%International Standard Atmosphere
if h<=11000
    T = 288.15 - 6.5*h*10^(-3);
    P = 101325*(1 - 0.2481*h/11000)^(5.2555);
    rho = P/(287*T);
else
    T = 216.65;
    P = 101325*0.2234*exp(-1.7344*(h/11000 - 1));
    rho = P/(287*T);
end
end

