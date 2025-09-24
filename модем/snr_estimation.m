% clear, clc

EbN0 = 13;
chan = comm.AWGNChannel('EbNo',EbN0,'BitsPerSymbol',1);
scale = 1;

for i = 1:1000
    bits = randi([0 1], 1280, 1);
    txSyms = pskmod(bits, 4, "InputType","bit");
    rxSyms = chan(txSyms);
    rxSyms = scale*rxSyms;
    
%     rxEb = (abs(mean(rxSyms(real(rxSyms*exp(1j*pi/4))>0 & imag(rxSyms*exp(1j*pi/4))>0))))^2;
%     rxVar = var(rxSyms(real(rxSyms*exp(1j*pi/4))>0 & imag(rxSyms*exp(1j*pi/4))>0));

%     rxEb = (abs(mean(rxSyms(real(rxSyms*exp(1j*pi/4))>0 & imag(rxSyms*exp(1j*pi/4))>0))))^2 + ...
%         (abs(mean(rxSyms(real(rxSyms*exp(1j*pi/4))<0 & imag(rxSyms*exp(1j*pi/4))>0))))^2 + ...
%         (abs(mean(rxSyms(real(rxSyms*exp(1j*pi/4))<0 & imag(rxSyms*exp(1j*pi/4))<0))))^2 + ...
%         (abs(mean(rxSyms(real(rxSyms*exp(1j*pi/4))>0 & imag(rxSyms*exp(1j*pi/4))>0))))^2;
%     rxVar = var(rxSyms(real(rxSyms*exp(1j*pi/4))>0 & imag(rxSyms*exp(1j*pi/4))>0)) + ...
%         var(rxSyms(real(rxSyms*exp(1j*pi/4))<0 & imag(rxSyms*exp(1j*pi/4))>0)) + ...
%         var(rxSyms(real(rxSyms*exp(1j*pi/4))<0 & imag(rxSyms*exp(1j*pi/4))<0)) + ...
%         var(rxSyms(real(rxSyms*exp(1j*pi/4))>0 & imag(rxSyms*exp(1j*pi/4))<0));

%     EbN0_est = 10*log10(rxEb/rxVar);
%     E2 = mean(abs(rxSyms).^2);
%     E4 = mean(abs(rxSyms).^4);
% 
%     a = 1;
%     b = -2*E2;
%     c = E4 - E2.^2;
%     D = b.^2 - 4*a*c;
% 
%     if D < 0 
%         continue
%     end
% 
%     Pn1 = (-b+sqrt(D))/2*a;
%     Pn2 = (-b-sqrt(D))/2*a;
% 
%     Ps1 = E2 - Pn1;
%     Ps2 = E2 - Pn2;
% 
%     if Ps1 > 0
%         Ps = Ps1;
%         Pn = Pn1;
%     else
%         Ps = Ps2;
%         Pn = Pn2;
%     end
% 
%     EbN0_est = 10*log10(Ps/Pn);
%     test(i) = EbN0_est;
    test(i) = snr_est(rxSyms);
end

disp(mean(test))

