clear, clc

N = 1024;
cpLen = 0.25*N;
Fs = 10e6;
mapperType = 'qam'; %'psk', 'qam'
M = 4;
log2M = log2(M);
chanType = 'eva'; % 'eva', 'simple'
eqType = 'zf'; %'zf', 'mmse', 'non' = AWGN

nul = 0;
mtrx = DFnTmtrx(N - nul);
imtrx = ctranspose(mtrx);

switch mapperType
    case 'psk'
        constellation = pskmod(0:M-1, M, 0);
    case 'qam'
        constellation = qammod(0:M-1, M);
end

Es = mean(abs(constellation).^2);
Eb = Es/log2M;

%==============Параметры моделирования=============
CFO = 0;
dF = CFO*(Fs/N);
EbN0 = 0:40;
isNoise = 1;
numFrames = 100;
isCP = 1;
if ~isCP
    cpLen = 0;
end

% phaseRotate = diag(N-1, N-1);
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


% всякие параметры частотно-избиратльного канала
switch chanType
    case 'eva'
        chan.PDelaysNs = [0, 30, 150, 310, 370, 710, 1090, 1730, 2510] * 1e-9;
        chan.APGdB = [0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12, -16.9];
        chan.MDSHz = 0;
        chan.Seed = 42;
    case 'simple'
%         chan.PDelaysNs = [0 600 1200 1800 2400 3000 3600 4200 4800 5400] * 1e-9;
%         chan.APGdB = [0 0 0 0 0 0 0 0 0 0];
        chan.PDelaysNs = [0 600 1200 1800 2400] * 1e-9;
        chan.APGdB = [0 0 0 0 0];
        chan.MDSHz = 0;
        chan.Seed = 42;
end

rayleighChan = comm.RayleighChannel("SampleRate", Fs,...
                    "PathDelays", chan.PDelaysNs,...
                    "AveragePathGains", chan.APGdB,...
                    "MaximumDopplerShift", chan.MDSHz,...
                    "ChannelFiltering",1,...
                    "PathGainsOutputPort",true,...
                    "Visualization","off", ...
                    'RandomStream', 'mt19937ar with seed',...
                    'Seed', chan.Seed);


%============== НАЧАЛО МОДЕЛИРОВАНИЯ ================
EbN0 = [80 EbN0];
ofdm_res = zeros(length(EbN0), 1);
ocdm_res = zeros(length(EbN0), 1);
tic; 
for SNR = 1:length(EbN0)
    for Fr=1:numFrames % цикл по кадрам
        if SNR == 1  
            preambule_ofdm = cell(1, 2);
            preambule_ocdm = cell(1, 2);
        end
    %=============Transmitter==============
        txBits = randi([0 1], log2M*(N - nul), 1);
        data = bit2int(txBits, log2M);
        switch mapperType
            case 'psk'
                txSyms = pskmod(data, M, 0, InputType='integer');
            case 'qam'
                txSyms = qammod(data, M, "gray");
        end
        % OCDM
        txSignal_ocdm = imtrx*txSyms;
        if nul ~= 0
            txSignal_ocdm = resample(txSignal_ocdm, 2, 1);
        end
        if SNR == 1 
            preambule_ocdm{1, 2} = txSignal_ocdm;
        end
    
        % OFDM
        txSyms_ofdm = [zeros(nul/2, 1); txSyms; zeros(nul/2, 1)];
        txSignal_ofdm = sqrt(N)*ifft(ifftshift(txSyms_ofdm), N);
        if SNR == 1 
            preambule_ofdm{1, 2} = txSyms_ofdm;
        end
    
        % добавление CP
        if isCP
            tmp1 = [txSignal_ocdm(end - cpLen+1:end); txSignal_ocdm];
            tmp2 = [txSignal_ofdm(end - cpLen+1:end); txSignal_ofdm];
        
            txSignal_ocdm = tmp1;
            txSignal_ofdm = tmp2;
        
            clear tmp1 tmp2 tmp3;
        end
    
%         for SNR = 1:length(EbN0)
    %=============CHANNEL=====================   
            N0 = (Eb/10^(EbN0(SNR)/10));
    
            if isNoise
                noise = ((1+1j)*randn(size(txSignal_ocdm)))*sqrt(N0/2); 
            else
                noise = zeros(size(txSignal_ocdm));
            end
            
            NSR = N0/Eb; 
            
            if eqType == "non"
                rxSignal_ocdm_before_noise = txSignal_ocdm;
                rxSignal_ofdm_before_noise = txSignal_ofdm;
            else
                CR = rayleighChan([1; zeros(49, 1)]); % определение отклика канала
                rayleighChan.release();
                [~, index] = max(abs(CR)); % позицию самого сильного луча
   
                rxSignal_ocdm_before_noise = rayleighChan([txSignal_ocdm; zeros(index,1)]);
                rxSignal_ocdm_before_noise = rxSignal_ocdm_before_noise(index:index+N+cpLen-1);
                rayleighChan.release();
                rxSignal_ofdm_before_noise = rayleighChan([txSignal_ofdm; zeros(index,1)]);
                rxSignal_ofdm_before_noise = rxSignal_ofdm_before_noise(index:index+N+cpLen-1);
                rayleighChan.release();
            end
    
            % Добавление Шума
            switch nul
                case 0
                    rxSignal_ocdm = (rxSignal_ocdm_before_noise + noise).*exp(1i*2*pi*dF*(0:length(txSignal_ocdm) - 1)/Fs)';
                otherwise
                    rxSignal_ocdm = (rxSignal_ocdm_before_noise + noise*2).*exp(1i*2*pi*dF*(0:length(txSignal_ocdm) - 1)/Fs)';
            end
            rxSignal_ofdm = (rxSignal_ofdm_before_noise + noise).*exp(1i*2*pi*dF*(0:length(txSignal_ofdm) - 1)/Fs)';
    %==================RECIEVER========================
    %   Удаление CP
            
            if isCP
                tmp1 = rxSignal_ocdm(cpLen+1:end);
                tmp2 = rxSignal_ofdm(cpLen+1:end);
    %             rxPreambule = rxPreambule(cpLen+1:end);
        
                rxSignal_ocdm = tmp1;
                rxSignal_ofdm = tmp2;
        
                clear tmp1 tmp2
            end

            if SNR == 1 
                preambule_ocdm{1, 1} = fft(rxSignal_ocdm);
                preambule_ocdm{1, 2} = fft(preambule_ocdm{1, 2});

                preambule_ofdm{1, 1} = rxSignal_ofdm;
                preambule_ofdm{1, 1} = (1/sqrt(N))*fftshift(fft(preambule_ofdm{1, 1}, N));

                chanEst_ofdm = preambule_ofdm{1, 1}./preambule_ofdm{1, 2};
                chanEst_ocdm = preambule_ocdm{1, 1}./preambule_ocdm{1, 2};
                continue
            end
            %================================//========================================
%     %-------OFDM-------
            rxSyms_ofdm = fftshift(fft(rxSignal_ofdm, N));
            rxSyms_ofdm = rxSyms_ofdm(nul/2+1:end-nul/2, 1);
    %-------Эквализация-
            switch eqType
                case 'zf'
                    rxSyms_ofdm = rxSyms_ofdm./chanEst_ofdm; % ZF
                case 'mmse'
                    rxSyms_ofdm = rxSyms_ofdm.*(conj(chanEst_ofdm)./(abs(chanEst_ofdm).^2 + NSR));   % MMSE
            end
            rxSyms_ofdm = (1/sqrt(N))*rxSyms_ofdm;
%                 figure(1);plot(rxSyms_ofdm, 'o'); xlim([-2 2]), ylim([-2 2]); pause(0.1)
%                 figure(2);plot(rxSyms_ocdm, 'o'); xlim([-2 2]), ylim([-2 2]); pause(0.1)
    %-------Демодулирование--
            switch mapperType
                case 'psk'
                    rxBits_ofdm = pskdemod(rxSyms_ofdm, M, 0, 'gray');
                case 'qam'
                    rxBits_ofdm = qamdemod(rxSyms_ofdm, M, 'gray');
            end
            rxBits_ofdm = int2bit(rxBits_ofdm, log2M);
    
    %-------OCDM-------
            if nul ~= 0
                rxSignal_ocdm = resample(rxSignal_ocdm, 1, 2);
            end
    
    %-------Эквализация------        
            rxSignal_ocdm_fft = fft(rxSignal_ocdm, N); % переход в частотную область
            rxSignal_ocdm_rot = rxSignal_ocdm_fft.*phaseRotate;
%             rxSignal_ocdm_rot = rxSignal_ocdm_fft;
            switch eqType
                case 'zf'
%                     rxSignal_ocdm_eq = rxSignal_ocdm_rot./fftshift(chanEst_ofdm); % оценка ofdm
                    rxSignal_ocdm_eq = rxSignal_ocdm_rot./chanEst_ocdm; % оценка ocdm
                case 'mmse'
                    rxSignal_ocdm_eq = rxSignal_ocdm_rot.*(conj(chanEst_ocdm)./abs(chanEst_ocdm).^2 + NSR);
                otherwise
                    rxSignal_ocdm_eq = rxSignal_ocdm_rot;
            end
%             rxSignal_ocdm_eq = rxSignal_ocdm_eq.*phaseRotate;
            rxSyms_ocdm = ifft(rxSignal_ocdm_eq, N);
%             figure(1); plot(rxSyms_ocdm, 'o');xlim([-2 2]); ylim([-2 2]); pause(0.5);
%             figure(1); stem(10*log10((angle(chanEst_ocdm)))); hold on
    %-------Демодулирование--
            switch mapperType
                case 'psk'
                    rxBits_ocdm = pskdemod(rxSyms_ocdm, M, 0, "gray");
                case 'qam'
                    rxBits_ocdm = qamdemod(rxSyms_ocdm, M, "gray");
            end
            rxBits_ocdm = int2bit(rxBits_ocdm, log2M);
            ocdm_res(SNR) = ocdm_res(SNR) + sum(rxBits_ocdm~=txBits)/length(txBits);

            ofdm_res(SNR) = ofdm_res(SNR) + sum(rxBits_ofdm~=txBits)/length(txBits);
    end
end
toc;

ocdm_res = ocdm_res/numFrames;
ofdm_res = ofdm_res/numFrames;

EbN0 = EbN0(2:end);
ocdm_res = ocdm_res(2:end);
ofdm_res = ofdm_res(2:end);
%% draw ber
figure()
semilogy(EbN0, ocdm_res, '+-')
hold on
semilogy(EbN0, ofdm_res, '-.')
hold on
switch mapperType
    case 'psk'
        bpsk_ber = berawgn(EbN0, 'psk', M, 'nondiff');
    case 'qam'
        bpsk_ber = berawgn(EbN0, 'qam', M);
end
semilogy(EbN0, bpsk_ber)
hold on
grid on
% xlim([1, 10]);
ylim([1e-7, 1])

theoryName = int2str(M) + "-" + upper(mapperType);
switch eqType 
    case "zf"
        legend('OCDM ZF', 'OFDM ZF', theoryName);
    case "mmse"
        legend('OCDM MMSE', 'OFDM MMSE', theoryName);
    otherwise
        legend('OCDM', 'OFDM', theoryName);
end