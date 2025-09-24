clc;
clear;

M = 4;
N = 128;
EbN0 = 5:1:10;
NumAver = 1e2;
% s = pskmod(randi([0 1], N, 1), M, 0,"gray","InputType","bit");
for snr = 1:length(EbN0)
    for i = 1:1:NumAver

        s = pskmod(randi([0 1], N, 1), M, 0,"InputType","bit");

        SNR = convertSNR(EbN0(snr),"ebno","snr","BitsPerSymbol", log2(M));

        r = awgn(s, SNR);
        r = randi([1 10], 1, 1)*r;

        E2 = mean(abs(r).^2);
        E4 = mean(abs(r).^4);

        a = 1;
        b = -2*E2;
        c = E4 - E2.^2;
        D = b.^2 - 4*a*c;

        if D < 0 
            continue
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

        m(i) = 10*log10(Ps/Pn);

    end
    m = nonzeros(m);
    estimate(snr) = mean(m);
end
figure
plot(estimate-3, linewidth = 2);
hold on;
plot(EbN0, linewidth = 2);
legend('Est','SNR', 'Location','northwest')