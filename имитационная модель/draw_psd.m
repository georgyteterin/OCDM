clear, clc

load('ocdm_for_psd.mat')
load('ofdm_for_psd.mat')

[psd_ocdm, freqs] = pwelch(resample(ocdmGlobal, 2, 1), [], [], [], 1e7, 'centered');
[psd_ofdm, freqs] = pwelch(resample(ofdmGlobal, 2, 1), [], [], [], 1e7, 'centered');

figure;
plot(freqs/1e6, 10*log10(psd_ocdm))
hold on
plot(freqs/1e6, 10*log10(psd_ofdm))
grid on

xl = xlabel('\it f \rm, MГц');
xl.FontName = 'Times New Roman';

yl = ylabel('Мощность/Частота, дБ/Гц');
yl.FontName = 'Times New Roman';

t = title('Спектральная плотность мощности');
t.FontName = 'Times New Roman';
t.FontSize = 14;

leg = legend('OCDM', 'OFDM');
leg.FontName = 'Times New Roman';
