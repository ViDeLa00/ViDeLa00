function [cineq, ceq] = nonlcon2BF(vars, fc, q) 

        
     OPR = fc.OPR;
     TIT = fc.TIT;
     BPR = vars(1);
     FPR = vars(2);
    
    T4t = TIT*273.15;
    P0t = fc.P0*(1+ (fc.g-1)*(fc.M0^2)/2)^(fc.g/(fc.g-1));
    if q.e02 ~= 1             %if q.e02 is not given, for M<1 the pressures are q.equal and for M>1 && M<5 the pressures are given by a relation from MIL-E-5007-D
        P2t = fc.P0*(1+ q.e02*(fc.g-1)*(fc.M0^2)/2)^(fc.g/(fc.g-1));
    elseif fc.M0>1 && fc.M0<5         %% only if supersonic
        P1t = P0t*(1-0.076*(fc.M0-1)^1.35);
        P2t = P1t*q.pi12;
    else
        P1t = P0t;
        P2t = P1t*q.pi12;
    end
    T2t = fc.T0*(1 + (fc.g-1)*(fc.M0^2)/2);
    Cp = fc.R*fc.g/(fc.g-1);
    Cpe = fc.Re*fc.ge/(fc.ge-1);
    
    f = (Cpe*T4t - Cp*T2t*(1 + (OPR^((fc.g-1)/fc.g) - 1)/q.e23))/(q.eq*fc.L - Cpe*TIT + Cp*T2t*(1 + (OPR^((fc.g-1)/fc.g) - 1)/q.e23));
    
    esp = (fc.g-1)/fc.g;
    espe = (fc.ge)/(fc.ge -1);
    
    P5t = P2t*OPR*((1 - Cp*T2t*(OPR^esp -1)/((1+f)*Cpe*q.e445*T4t*q.e23))^espe)*(1 - Cp*T2t*BPR*(FPR^esp - 1)/((1+f)*Cpe*q.e455*q.efan*(T4t - Cp*T2t*(OPR^esp - 1)/((1+f)*Cpe*q.e23))))^espe;
    T5t = T4t - T2t*Cp*(OPR^esp -1)/((1+f)*Cpe*q.e23) - T2t*BPR*Cp*(FPR^esp -1)/((1+f)*Cpe*q.efan);
    
    %Potrebbe succedere che l'espansione sia tale da avere Temperature di
    %uscita negative, quindi effettuiamo un check rispetto alla Temperatura
    %minima raggiungibile, che è quella del Turbojet dalla stazione "45"
    
    adapted = false;
    if P5t/fc.P0 <= ((fc.ge+1)/2)^espe
        adapted = true;
        Ps = fc.P0;
    end
    if adapted
            Ts = T5t*(fc.P0/P5t)^((fc.ge-1)/fc.ge);
    else
            Ts = T5t*2/(fc.ge+1);     
    end
    
    Ttj = (fc.P0^((fc.ge-1)/fc.ge))*(T4t - Cp*T2t*(OPR^esp -1)/((1+f)*Cpe*q.e23))/(((OPR*P2t)^((fc.ge-1)/fc.ge))*(1 - (Cp*T2t)*(OPR^esp -1)/((1+f)*q.e445*T4t*Cpe*q.e23)));
    
    %controllo se la pressione P5t è maggiore della pressione ambiente
    
    cineq = [fc.P0 - P5t; Ttj - Ts];
    ceq = [];
end    



