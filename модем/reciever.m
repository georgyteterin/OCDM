clear, clc
load("txData_2905_1523.mat")
% load("txData.mat")
N = 64;
fs = 1e7;
mtrx = DFnTmtrx(N); %нужно добавить сохранение всех параметров в .mat файл

mapper = comm.PSKModulator(4, 0, "BitInput",true);
demapper = comm.PSKDemodulator(4, 0, "BitOutput",true);

phaseRotate = zeros(N, 1);
for k=0:N-1
    switch mod(N, 2)
        case 0
            phaseRotate(k+1, 1) = exp(-1j*(pi/N)*k^2);
        otherwise
            phaseRotate(k+1, 1) = exp(-1j*(pi/N)*k*(k-1));
    end
end
clear k;

%============Параметры============
eq = 0;
sdr_type = "usrp"; %'usrp'/'hackrf'

% ============Чтение============
if sdr_type == "hackrf"
    rx = audioread("14-13-39_2199000000Hz.wav");
    rx_baseband = rx(:,1) + 1i*rx(:,2);
elseif sdr_type == "usrp"
    fd = fopen('C:\records_usrp\30_05\rec1.bin','rb');
%     fd = fopen('C:\Users\гоша\Documents\VKR\tx_record.bin','rb');
    y = fread(fd, 'int16');
    fclose(fd);
    rx_baseband = y(1:2:end) + 1i*y(2:2:end);
end

% rec1 ... delta f = +90
% rec2 ... delta f = + 123
% rec3 ... delta f = 132
% rec4 ... delta f = 149
% rec5 ... delta f = 90
% rx_baseband = rx_baseband(1:1e6);
% rx_baseband = rx_baseband.*exp(1i*2*pi*-1e6*(0:length(rx_baseband)-1)/fs).';
rx_baseband = rx_baseband.*exp(1i*2*pi*90*(0:length(rx_baseband)-1)/fs).';
rx_baseband = decimate(rx_baseband, 2, 'fir');
fs = fs/2;

corr = xcorr(rx_baseband, preambule1);
corr = corr(length(rx_baseband)-length(preambule1):end);
%%
gate = max(abs(corr))*5/6;
[val, pos] = findpeaks(abs(corr), "MinPeakHeight", gate);
%%
globalSyms = [];
EbN0_est = [];
SER = [];
d_phi_global = [];
phi_global = [];
chan_est_global = [];
for ind=1:length(pos)-1
    phi = angle(corr(pos(ind)));
    phi_global = [phi_global phi];
       
    shift = 0;

    rx_ocdm_frame = rx_baseband(pos(ind) + 2*length(preambule2)+shift:pos(ind) + 2*length(preambule2) + 880 - 1+shift);
    rx_ocdm_frame = rx_ocdm_frame*exp(1i*-phi);

    rx_pr2 = rx_baseband(pos(ind):pos(ind)+length(preambule2) - 1);
    rx_pr3 = rx_baseband(pos(ind)+length(preambule2):pos(ind)+2*length(preambule2)-1);
    rx_pr2 = rx_pr2*exp(1i*-phi);
    rx_pr3 = rx_pr3*exp(1i*-phi);
    cor_val2 = sum((rx_pr2).*conj(preambule2));
    cor_val3 = sum((rx_pr3).*conj(preambule2));

    d_phi = angle(cor_val3) - angle(cor_val2);
    d_phi_global = [d_phi_global d_phi];

    freq_shift(ind) = d_phi/(length(preambule2)/fs)/2/pi;
    rx_ocdm_frame = (2.5/max(abs(rx_ocdm_frame)))*rx_ocdm_frame;
%      pwelch(resample(rx_ocdm_frame, 2, 1), [], [], [], 2*fs, 'centered'), pause(0.01);
    rx_ocdm_frame = rx_ocdm_frame.*exp(1i*-freq_shift(ind)*(0:length(rx_ocdm_frame)-1)/fs).';
    % выделение пилотного сигнала
    pilot_ocdm_sym = rx_ocdm_frame(1:80);
    pilot_ocdm_sym = pilot_ocdm_sym(17:end);
    ref_pilot_no_CP = ref_pilot(17:end);
    chanEst = (fft(pilot_ocdm_sym, N).*phaseRotate)./(fft(ref_pilot_no_CP, N).*phaseRotate);
    chan_est_global = [chan_est_global chanEst];
    pilot_ocdm_sym_fft = fft(pilot_ocdm_sym, N);
    pilot_ocdm_sym_rot = pilot_ocdm_sym_fft.*phaseRotate;
    pilot_ocdm_sym_eq = pilot_ocdm_sym_rot./chanEst;
    pilot_syms = ifft(pilot_ocdm_sym_eq, N);

    %     chanEst = (fft(ref_pilot, N))./(fft(pilot_ocdm_sym, N));
    %     plot(mtrx*pilot_ocdm_sym, 'o'); pause(0.05);

    %     globalSyms = [];
    %     SymbolError = zeros(10, 1);
%     plot(10*log10(abs(fftshift(chanEst)))); pause(0.1);
%     pwelch(resample(pilot_ocdm_sym, 2, 1), [], [], [], fs*2, 'centered'); pause(0.1);
    for symNum=2:11% ДЛЯ ОЦЕНКИ КАНАЛА ПОМЕНЯТЬ НА 2:11
        symLen = 80;
        ocdm_sym = rx_ocdm_frame(1+symLen*(symNum-1):symLen*symNum);
%         pwelch(ocdm_sym, [], [], [], 10e6, 'centered'); pause(0.05)
        ocdm_sym = ocdm_sym(17:end);
%         ocdm_sym = ocdm_sym.*exp(1i*2*pi*-freq_shift*(0:length(ocdm_sym)-1)/fs).';
        ocdm_sym_fft = fft(ocdm_sym, N);
        ocdm_sym_rot = ocdm_sym_fft.*phaseRotate;
        if ~eq
            ocdm_sym_eq = ocdm_sym_rot;
        else
            ocdm_sym_eq = ocdm_sym_rot./chanEst;
        end
        syms = ifft(ocdm_sym_eq, N);
        syms = syms - mean(syms);
        %         plot(syms, 'o');xlim([-0.2 0.2]); ylim([-0.2 0.2]); pause(0.05);
        %         syms = syms*exp(-1j*pi/8);
        globalSyms = [globalSyms; syms];
        %         syms = mtrx*ocdm_sym*exp(1j*pi/2);
        

        EbN0_est = [EbN0_est snr_est(syms)]; 
        
%         if (snr_est(syms) > 5.5 && snr_est(syms) < 6.5)
%             globalSyms = [globalSyms; syms];
%         end

        bits = demapper(syms);
        BitError = biterr(bits, all_data(symNum-1, :)');
        SER = [SER BitError];
%         plot(syms, 'o'); xlim([-1.5 1.5]); ylim([-1.5 1.5]);pause(0.05);
    end
end

