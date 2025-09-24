load("results/awgn_channel/AWGN CFO=0 psk-4.mat")
theory = berawgn(EbN0, 'psk', 4, 'nondiff');
ocdm_res = theory;
ofdm_res = theory;
figure()
semilogy(EbN0, ocdm_res, 'b-d')
hold on
semilogy(EbN0, ofdm_res, 'r-s')
hold on
semilogy(EbN0, theory)
hold on
grid on
ylim([1e-7, 1])


xl = xlabel('Eb/N0, дБ');
xl.FontName = 'Times New Roman';
xl.FontSize = 12;
yl = ylabel('Вероятность битовой ошибки');
yl.FontName = 'Times New Roman';
yl.FontSize = 12;
l = legend('OCDM', 'OFDM', 'QPSK') ;
l.FontName = 'Times New Roman';
l.FontSize = 10;

%%
load("results/awgn_channel/y.mat")
theory = berawgn(EbN0, 'psk', 4, 'nondiff');
figure()
semilogy(EbN0, ocdm_res, 'b-d')
clear ocdm_res
hold on
load("results/awgn_channel/n.mat")
semilogy(EbN0, ocdm_res, 'r-s')
hold on
semilogy(EbN0, theory, 'g-s')
hold on
grid on
ylim([1e-7, 1])


xlabel('Eb/N0')
ylabel('BER')
legend('upsampled', 'clear', 'theory')