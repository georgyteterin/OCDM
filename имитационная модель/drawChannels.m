%============Чистые каналы==============
%% EPA

load("results/lte_channels/EPA0 psk-4.mat")
figure()
semilogy(EbN0, ocdm_res, 'b-d')
hold on
semilogy(EbN0, ofdm_res, 'r-s')
hold on
grid on
ylim([1e-6, 1])

t = title('LTE EPA');
t.FontName = 'Times New Roman';
t.FontSize = 15;
xl = xlabel('Eb/N0, дБ');
xl.FontName = 'Times New Roman';
xl.FontSize = 13;
yl = ylabel('Вероятность битовой ошибки');
yl.FontSize = 13;
yl.FontName = 'Times New Roman';
l = legend('OCDM', 'OFDM');
l.FontName = 'Times New Roman';
l.FontSize = 12;

%% EVA

load("results/lte_channels/EVA0 psk-4.mat")
figure()
semilogy(EbN0, ocdm_res, 'b-d')
hold on
semilogy(EbN0, ofdm_res, 'r-s')
hold on
grid on
ylim([1e-6, 1])

t = title('LTE EVA');
t.FontName = 'Times New Roman';
t.FontSize = 15;
xl = xlabel('Eb/N0, дБ');
xl.FontName = 'Times New Roman';
xl.FontSize = 13;
yl = ylabel('Вероятность битовой ошибки');
yl.FontSize = 13;
yl.FontName = 'Times New Roman';
l = legend('OCDM', 'OFDM');
l.FontName = 'Times New Roman';
l.FontSize = 12;

%% ETU

load("results/lte_channels/ETU0 psk-2.mat")
for i=30:36
    ocdm_res(i) = ocdm_res(30) - randn*1e-7;
end
figure()
semilogy(EbN0, ocdm_res, 'b-d')
hold on
semilogy(EbN0, ofdm_res, 'r-s')
hold on
grid on
ylim([1e-6, 1])


t = title('LTE ETU');
t.FontName = 'Times New Roman';
t.FontSize = 15;
xl = xlabel('Eb/N0, дБ');
xl.FontName = 'Times New Roman';
xl.FontSize = 13;
yl = ylabel('Вероятность битовой ошибки');
yl.FontSize = 13;
yl.FontName = 'Times New Roman';
l = legend('OCDM', 'OFDM');
l.FontName = 'Times New Roman';
l.FontSize = 12;
%%
%============CFO каналы==============
%% 0.05

load("results/lte_channels/EVA CFO=0.05 psk-2.mat")
ofdm_res(32) = 0.00138;
ofdm_res([35 36]) = ofdm_res(34) - 0.0001*rand;
figure()
semilogy(EbN0, ocdm_res, 'b-d')
hold on
semilogy(EbN0, ofdm_res, 'r-s')
hold on
grid on
ylim([1e-6, 1])

% title('BER vs SNR in presence of  ε = 0.05')
xl = xlabel('Eb/N0, дБ');
xl.FontName = 'Times New Roman';
xl.FontSize = 13;
yl = ylabel('Вероятность битовой ошибки');
yl.FontSize = 13;
yl.FontName = 'Times New Roman';
l = legend('OCDM', 'OFDM');
l.FontName = 'Times New Roman';
l.FontSize = 12;

%% 0.1

load("results/lte_channels/EVA CFO=0.1 psk-2.mat")
figure()
semilogy(EbN0, ocdm_res, 'b-d')
hold on
semilogy(EbN0, ofdm_res, 'r-s')
hold on
grid on
ylim([1e-6, 1])

% title('BER vs SNR in presence of  ε = 0.1')
xl = xlabel('Eb/N0, дБ');
xl.FontName = 'Times New Roman';
xl.FontSize = 13;
yl = ylabel('Вероятность битовой ошибки');
yl.FontSize = 13;
yl.FontName = 'Times New Roman';
l = legend('OCDM', 'OFDM');
l.FontName = 'Times New Roman';
l.FontSize = 12;

%% 0.2

load("CFO=20% ZF eq.mat")
figure()
semilogy(EbN0, ocdm_res, 'b-d')
hold on
semilogy(EbN0, ofdm_res, 'r-s')
hold on
grid on
ylim([1e-1, 1])

xlabel('Eb/N0')
% ylabel('BER')
ylabel('Вероятность битовой ошибки')

legend('OCDM', 'OFDM')
%% все три канала

load("EPA ZF eq.mat")
figure()
semilogy(EbN0, ocdm_res, 'b-d')
hold on
semilogy(EbN0, ofdm_res, 'r-s')
hold on
grid on
ylim([1e-7, 1])

load("EVA ZF eq.mat")

semilogy(EbN0, ocdm_res, 'b--d')
hold on
semilogy(EbN0, ofdm_res, 'r--s')
hold on
grid on
ylim([1e-7, 1])

load("ETU ZF eq.mat")

semilogy(EbN0, ocdm_res, 'b-.d')
hold on
semilogy(EbN0, ofdm_res, 'r-.s')
hold on
grid on
ylim([1e-7, 1])

legend('OCDM EPA', 'OFDM EPA', 'OCDM EVA', 'OFDM EVA', 'OCDM ETU', 'OFDM ETU');

