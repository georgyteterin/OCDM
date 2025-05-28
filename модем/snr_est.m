function snr = snr_est(r)
    E2 = mean(abs(r).^2);
    E4 = mean(abs(r).^4);
    
    a = 1;
    b = -2*E2;
    c = E4 - E2.^2;
    D = b.^2 - 4*a*c;
    
    if D < 0 
        error("D<0")
%         D = 1;
    end
    
    Pn1 = (-b+sqrt(D))/2*a;
    Pn2 = (-b-sqrt(D))/2*a;
    
    Ps1 = E2 - Pn1;
    Ps2 = E2 - Pn2;
    
    if Ps1 > 0
        Ps = Ps1;
        Pn = Pn1;
    else
        Ps = Ps2;
        Pn = Pn2;
    end
    
    snr = 10*log10(Ps/Pn);
end