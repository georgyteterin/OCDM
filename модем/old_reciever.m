clear, clc
load("txData.mat")
ref_pilot = pilot_ocdm_symbol;
clear pilot_ocdm_symbol
N = 64;
mtrx = DFnTmtrx(N); %нужно добавить сохранение всех параметров в .mat файл
% mapper = comm.RectangularQAMModulator(4, "BitInput",true);
mapper = comm.PSKModulator(4, 0, "BitInput",true);
% demapper = comm.RectangularQAMDemodulator(4, "BitOutput",true);
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
eq = 1;

sdr_type = "usrp"; %'usrp'/'hackrf'

% load("rx_baseband.mat")

if sdr_type == "hackrf"
    rx = audioread("14-13-39_2199000000Hz.wav");
    rx_baseband = rx(:,1) + 1i*rx(:,2);
elseif sdr_type == "usrp"
    fd = fopen('C:\records_usrp\28_05\rec2.bin','rb');
%     fd = fopen('C:\Users\гоша\Documents\VKR\tx_record.bin','rb');
    y = fread(fd, 'int16');
    fclose(fd);
    rx_baseband = y(1:2:end) + 1i*y(2:2:end);
end
rx_baseband = rx_baseband(100:end);
% rx_baseband = rx_baseband - mean(rx_baseband);
[spect_values, spect_pos] = pwelch(rx_baseband.^2, [],[],[],1e7, 'centered');
[max_spec, max_pos] = max(spect_values);
freq_shift = spect_pos(max_pos);
rx_baseband = rx_baseband.*exp(1i*2*pi*-freq_shift/2*(0:length(rx_baseband)-1)/1e7).';


% rx_dec = decimate(rx_baseband, 2, 'fir');
corr = xcorr(rx_baseband, preambule);
corr = corr(length(rx_baseband)-length(preambule):end);
%%
gate = max(abs(corr))*3/4;
[val, pos] = findpeaks(abs(corr), "MinPeakHeight", gate);

%%
globalSyms = [];
EbN0_est = [];
EbN0_paket = [];
SER = [];
massive_phi = [];
carrSynch = comm.CarrierSynchronizer('Modulation','QPSK', 'SamplesPerSymbol', 1, 'ModulationPhaseOffset', 'Custom', 'CustomPhaseOffset', 0);
for ind=1:length(pos)
    num = ind;
    phi = angle(corr(pos(num)));
%     massive_phi = [massive_phi phi];
    shift = 0;
    rx_ocdm_frame = rx_baseband(pos(num)+shift:pos(num)+10001*80*2+shift-1);
    rx_ocdm_frame = rx_ocdm_frame*exp(1i*-(angle(corr(pos(num)))));
    pwelch(rx_ocdm_frame, [], [], [], 10e6, 'centered'); pause(0.05)
%     rx_ocdm_frame = 2.5*rx_ocdm_frame/max(abs(rx_ocdm_frame));
%     pwelch(rx_ocdm_frame, [],[],[],2e7, 'centered'); pause(0.1);
    tmp = movmean(rx_ocdm_frame, 300);
    tmp = rx_ocdm_frame - tmp;
    tmp = rx_ocdm_frame;
    %     rx_ocdm_frame_dec = decimate(rx_ocdm_frame, 2, 'fir');
    tmp_dec = decimate(tmp, 2, 'fir');
    

    % выделение пилотного сигнала
    pilot_ocdm_sym = tmp_dec(1:80);
    pilot_ocdm_sym = pilot_ocdm_sym(17:end);
    ref_pilot_no_CP = ref_pilot(17:end);
    chanEst = (fft(pilot_ocdm_sym, N).*phaseRotate)./(fft(ref_pilot_no_CP, N).*phaseRotate);
%     chanEst(1) = 1;
    pilot_ocdm_sym_fft = fft(pilot_ocdm_sym, N);
    pilot_ocdm_sym_rot = pilot_ocdm_sym_fft.*phaseRotate;
    pilot_ocdm_sym_eq = pilot_ocdm_sym_rot./chanEst;
    pilot_syms = ifft(pilot_ocdm_sym_eq, N);
    %     chanEst = (fft(ref_pilot, N))./(fft(pilot_ocdm_sym, N));
    %     plot(mtrx*pilot_ocdm_sym, 'o'); pause(0.05);

    %     globalSyms = [];
    %     SymbolError = zeros(10, 1);
%     plot(10*log10(abs(fftshift(chanEst)))); hold on; pause(0.1);
    for symNum=7000:10001 % ДЛЯ ОЦЕНКИ КАНАЛА ПОМЕНЯТЬ НА 2:11
        symLen = 80;
        ocdm_sym = tmp_dec(1+symLen*(symNum-1):symLen*symNum);
        pwelch(ocdm_sym, [], [], [], 10e6, 'centered'); pause(0.05)
        ocdm_sym = ocdm_sym(17:end);
        phi_sym = angle(ocdm_sym(1));
        massive_phi = [massive_phi phi_sym];
        ocdm_sym_fft = fft(ocdm_sym, N);
        ocdm_sym_rot = ocdm_sym_fft.*phaseRotate;
        if ~eq
            ocdm_sym_eq = ocdm_sym_rot;
        else
            ocdm_sym_eq = ocdm_sym_rot./chanEst;
        end
        syms = ifft(ocdm_sym_eq, N);
        syms = syms - mean(syms);
        syms = carrSynch(syms);
        %         plot(syms, 'o');xlim([-0.2 0.2]); ylim([-0.2 0.2]); pause(0.05);
        %         syms = syms*exp(-1j*pi/8);
        globalSyms = [globalSyms; syms];
        %         syms = mtrx*ocdm_sym*exp(1j*pi/2);
        rxEb = (abs(mean(syms(real(syms*exp(1j*pi/4))>0 & imag(syms*exp(1j*pi/4))>0))))^2;
        rxVar = var(syms(real(syms*exp(1j*pi/4))>0 & imag(syms*exp(1j*pi/4))>0));
        

        EbN0_est = [EbN0_est snr_est(syms)-3]; 


        bits = demapper(syms);
        BitError = biterr(bits, all_data(symNum-1, :)');
        SER = [SER BitError];
%         plot(syms, 'o'); xlim([-2 2]); ylim([-2 2]);pause(0.05);
    end
    EbN0_paket = [EbN0_paket mean(EbN0_est((ind-1)*10+1:ind*10))];
%     carrSynch.release;
%     carrSynch.reset;
    %     meanSymbolError = mean(SymbolError);
    %     frameError(ind, 1) = meanSymbolError;
end
% meanFrameError = mean(frameError);
% close all


