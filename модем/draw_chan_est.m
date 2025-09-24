clear, clc

load chan_est_rec3.mat 

abs(mean(chan_est_global(:, 1:50), 2));
% plot(10*log10(fftshift(abs(ans))))
plot(fftshift(ans));
grid on
xl = xlabel('отсчеты бпф');
yl = ylabel('коэффициент передачи');
xl.FontName = 'Times New Roman';
yl.FontName = 'Times New Roman';

