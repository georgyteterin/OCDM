
clc; clear; close all

% Параметры
    params.NFFT = 2^8;                     % размер сообщения/NFFT
    params.NumSubcar = 52;               % Число поднесущих "пилот + данные"
    params.NumData = 48;
    params.NumPilots = 4;
    params.T = 4e-6;                        % длительность символа включая префикс в с
    params.cpDur = 0.8e-6;                 % длительность префикса в с    
    params.M = 16;                          % Размер созвездия для информационных поднесущих
    h2dB = 1:100;
    BER = zeros(size(h2dB));

% Задания
    eq = 0;
    cPref = 1;
    cCycle = 1;
    rayCh = 0;
    dft = 1;
% Вычисляемые параметры
    log2M = log2(params.M);

    params.Fs = (5/4)*params.NFFT/params.T;       % частота дискретизации
    params.cpLen = params.cpDur * params.Fs;
    params.left = (params.NFFT - params.NumSubcar) / 2;       % Размер защитного интервала слева
    params.right = (params.NFFT - params.NumSubcar) / 2;      % Размер защитного интервала справа
    params.infBitsNum = (params.NumData - 1)*log2M;

%% Передающая часть 
% Вычисление индексов информац-х, пилотных и защитных поднесущих
    % Индексы защитных интервалов
    ind.left = (1 : params.left).';
    ind.right = ((params.NFFT - params.right + 1) : params.NFFT).';

    % Индексы пилотных поднесущих
    ind.pilots = [1 14 28 40] + ind.left(end);

    % Индекс главной поднесущей
    ind.mainCarrier = params.NFFT/2 + 1;
    
    % Индексы информационных поднесущих
    ind.data = (ind.left(end)+1 : (ind.right(1)-1)).';
    clear pilotInds1 pilotInds2;

    % Исключение индексов поднесущих, которые являются пилотными или главной.
    tmplogic = false(size(ind.data));

    for i = 1 : length(ind.data)
        for j = 1 : length(ind.pilots)
            if (ind.data(i) == ind.pilots(j)) || (ind.data(i) == ind.mainCarrier)
                tmplogic(i) = 1;
                break
            end
        end
    end
    ind.data(tmplogic) = [];
    clear tmplogic i j;
    

for SNR = 1:length(h2dB)
    for Frame = 1:20
    % Формирование OFDM сигнала
    % Генерация информационной последовательности
        txBits = randi([0 1], params.infBitsNum, 1);
    
    % При необходимости дополнение информационной последовательности нулями
        if mod(length(txBits), (params.NumData - 1)* log2M) ~= 0
            txBitsExtended = [txBits;...
                               zeros(params.data * log2M - mod(length(txBits), params.data * log2M), 1)];
        else
            txBitsExtended = txBits;
        end
    
    % Формировнаие мод-нных символов информационных поднесущих
        tmp = reshape(txBitsExtended, (params.NumData - 1)*log2M, []);
        txSymbols = qammod(tmp, params.M, 'gray', 'InputType', 'bit');
        constellation = unique(txSymbols);
        Ps = mean(abs(constellation).^2);
        Pb = Ps/log2M;
        clear tmp;
    % Формировнаие мод-нных символов пилотных поднесущих
    %     pilotSymbs = qammod(mod(0:params.pilot-1, MPilot).', MPilot, 'gray', 'UnitAveragePower', true, 'InputType', 'integer');
        pilotSymbs = 1+1j;
    
    % Заготовка формирования OFDM сигнала
        txMatrix = zeros(params.NFFT, size(txSymbols, 2));
    
    % Заполнение заготовки символами
        % Цикл по символам
        for i = 1 : size(txMatrix, 2)
            % Заполнение пилотами
            txMatrix(ind.pilots, i) = pilotSymbs;
    
            % Заполнение информацией
            txMatrix(ind.data, i) = txSymbols(:, i);
        end
    
    % Обратное дискретное преобразование Фурье
 %         txSignal = txMatrix;
     if dft
        txSignal = ifft(ifftshift(txMatrix), params.NFFT);
     else
        txSignal = txMatrix;
     end
    %     txSignal = IDFnT(txMatrix);
        
    % Защитный интервал 
     if(cPref == 1)
         txSignalWithCp = zeros(params.cpLen + params.NFFT, size(txSignal, 2));
         if(cCycle == 1)
            txSignalWithCp(1 : params.cpLen, :) = txSignal(end - params.cpLen + 1 : end, :);
         end
         txSignalWithCp(params.cpLen + 1 : end, :) = txSignal;
    
    % Передаваемый сигнал
        txSignal = txSignalWithCp(:);
    else
        txSignal = txSignal(:);
    end
    clear txSignalWithCp;
    
    % Канал связи
        % Релеевский канал (EVA)
        if(rayCh == 1)
            EVA.PDelaysNs = [0, 30, 150, 310, 370, 710, 1090, 1730, 2510] * 1e-9;
            EVA.APGdB = [0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12, -16.9];
            EVA.MDSHz = 0;
    %         rayleighChan = comm.RayleighChannel("SampleRate", Fs, "PathDelays", EVA.PDelaysNs, "AveragePathGains", EVA.APGdB, "MaximumDopplerShift", EVA.MDSHz,'Visualization','Impulse and frequency responses');
            rayleighChan = comm.RayleighChannel("SampleRate", params.Fs, "PathDelays", EVA.PDelaysNs, "AveragePathGains", EVA.APGdB, "MaximumDopplerShift", EVA.MDSHz);
            rxSignal = rayleighChan(txSignal);
            rayleighChan.release();
            rxSignalBeforeAWGN = rxSignal;
        else
            rxSignalBeforeAWGN = txSignal;
        end
        % АБГШ
        Sigma = sqrt(Pb*10^(-SNR/10) / 2);
        noise = randn(length(rxSignalBeforeAWGN), 2) * [1; 1i] * Sigma;
%         noise = zeros(size(rxSignalBeforeAWGN));
        rxSignal = rxSignalBeforeAWGN + noise;
      
    % Приёмная часть
    % Обработка принятого сигнала
        % Разбиение по кадрам
        if (cPref == 0)
            rxMatrix = reshape(rxSignal, params.NFFT, []);
        else
            rxMatrix = reshape(rxSignal, params.NFFT + params.cpLen, []);
            rxMatrix = rxMatrix(params.cpLen + 1 : end, :);
        end
        % Прямое ДПФ
        if dft
            rxMatrix = fftshift(fft(rxMatrix,  params.NFFT));
        end
    %     rxMatrix = DFnT(rxMatrix);
    
    % Извлечение символов
        % Заготовка под пилоты
        rxPilots = zeros(params.NumPilots, size(rxMatrix, 2));
    
        % Заготовка под извлечённую информацию
        rxSymbols = zeros(params.NumData - 1, size(rxMatrix, 2));
    
        % Цикл по символам 
        for i = 1 : size(rxMatrix, 2)
            % Извлечение мод-нных символов
                rxSymbols(:, i) = rxMatrix(ind.data, i);
    
            % Извлечение пилотов
            rxPilots(:, i) = rxMatrix(ind.pilots, i);
        end
    
    % Эквалайзер
    if(eq == 1)
        rxPilots1D = rxPilots(:, 1);
        txPilots1D = pilotSymbs(:);
        H =  rxPilots1D ./ txPilots1D;
        Heq = interp1(ind.pilots, H, ind.data, 'spline');
        rxSymbols = rxSymbols ./ Heq;
    end
    % Демодуляция
        rxBitsExtended = qamdemod(rxSymbols, params.M, 'gray', 'OutputType', 'bit');
        rxBitsExtended = rxBitsExtended(:);
    
        rxBits= rxBitsExtended(1:length(txBits));
    
    % Вычисление BER
        BER(SNR) = BER(SNR) + sum((txBits ~= rxBits));
    end
end

BER = BER/(length(txBits)*Frame);
%% draw BER
BERth = berawgn(h2dB, "qam", params.M);

% berfit(h2dB, BER);
semilogy(BER);
grid on;
hold on;
semilogy(BERth);
ylim([10^-6 1]);
% title("BER comparison");
