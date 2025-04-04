clc; clear; 
% close all


% Параметры
    params.N = 1024;                        % размер сообщения/NFFT
    params.NumSubcar = 1024;                % Число поднесущих "пилот + данные"
    params.pilotStep = 3;
    params.T = 3e-6;                        % длительность символа включая префикс в с
    params.Fs = 10e6;   
    params.M = 16;                          % Размер созвездия для информационных поднесущих
    params.Mp = 4;                          % Размер созвездия для пилотных поднесущих
    params.FreqShift = 1/1e5;
    h2dB = 0:40;
    BER = zeros(size(h2dB));
    NumFrame = 400;
    



    tmp_sig = [];
% Задания
    eq = 0;
    prefix = 'CP'; % 'CP'/'ZP'
    rayCh = 0;
    params.style = ["ocdm", "ofdm"];

% Вычисляемые параметры
    log2M = log2(params.M);
    log2Mp = log2(params.Mp);
    dM = log2M/log2Mp;
    
    params.NumPilots = floor(params.NumSubcar/(params.pilotStep-1));
    params.cpLen = 0.25*params.N;
    params.left = (params.N - params.NumSubcar) / 2;       % Размер защитного интервала слева
    params.right = (params.N - params.NumSubcar) / 2;      % Размер защитного интервала справа

    params.DFnT_matrix = DFnTmtrx(params.N);
    params.IDFnT_matrix = ctranspose(params.DFnT_matrix);

    if ~rayCh
        params.ChanType = 'awgn';
    else
        params.ChanType = 'eva';
    end
for st = 1:length(params.style)
    %% параметры сохранения
        params.SaveFileName = num2str(params.style(st)) + " freq shift = " + num2str( ...
            params.FreqShift) + "pi";
    %     params.SaveDirName = "";
        if ~isfield (params, 'SaveDirName')
            params.SaveDirName = "DefaultResults";
        end
        if ~isfolder(params.SaveDirName)
            mkdir(params.SaveDirName);
        end
        params.SaveDirName = params.SaveDirName + '/';
        if ~isfield (params, 'SaveFileName')
            params.SaveFileName = params.style(st) + ' '+ params.ChanType + '.mat';
        else
            params.SaveFileName = params.SaveDirName + params.SaveFileName + '.mat';
        end
    %% Передающая часть 
    % Вычисление индексов информац-х, пилотных и защитных поднесущих
        % Индексы защитных интервалов
        if params.NumSubcar == params.N
            ind.left = 0;
            ind.right = params.N + 1;
        else
            ind.left = (1 : params.left).';
            ind.right = ((params.N - params.right + 1) : params.N).';
        end
        
        % Индекс главной поднесущей
        ind.mainCarrier = params.N/2 + 1;
    
        % Индексы пилотных поднесущих
        ind.pilots(1) = 1;
        for i=2:params.NumPilots-1
            ind.pilots(i) = ind.pilots(i-1)+params.pilotStep;
            if ind.pilots(i) >= ind.right
                ind.pilots(i) = [];
                break;
            end
        end
        params.NumPilots = length(ind.pilots);
        params.NumData = params.NumSubcar - params.NumPilots;
        params.infBitsNum = (params.NumData - 1)*log2M;
        ind.pilots(params.NumPilots) = params.NumSubcar;
        ind.pilots = ind.pilots + ind.left(end);
    
          
        % Индексы информационных поднесущих
        ind.data = (ind.left(end)+1 : (ind.right(1)-1)).';
        clear pilotInds1 pilotInds2;
    
        % Исключение индексов поднесущих, которые являются пилотными или главной.
        tmplogic = false(size(ind.data));
    
        for i = 1 : length(ind.data)
            if ind.data(i) == ind.mainCarrier
                tmplogic(i) = 1;
            end
            for j = 1 : length(ind.pilots)
                if (ind.data(i) == ind.pilots(j)) || (ind.data(i) == ind.mainCarrier)
                    tmplogic(i) = 1;
                    break
                end
            end
        end
        ind.data(tmplogic) = [];
        clear tmplogic i j;
        
    
    for s = 1:length(h2dB)
        SNR = h2dB(s);
        EbN0 = convertSNR(SNR,'ebno',BitsPerSymbol=log2M);
        for Frame = 1:NumFrame
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
            txSymbols = qammod(txBitsExtended, params.M, 'gray', 'UnitAveragePower', true, 'InputType', 'bit');
            constInt = 0:params.M - 1;
            constellation = qammod(constInt, params.M);
    
            Eb = 1/params.N/log2M;
            clear tmp;
        % Формировнаие мод-нных символов пилотных поднесущих
            txPilots = qammod(mod(0:params.NumPilots-1, 4).', params.Mp, 'gray', 'UnitAveragePower', true, 'InputType', 'integer');
    
        % Заготовка формирования OFDM сигнала
            txMatrix = zeros(params.N, size(txSymbols, 2));
        
        % Заполнение заготовки символами
            % Цикл по символам
            for i = 1 : size(txMatrix, 2)
                % Заполнение пилотами
                txMatrix(ind.pilots, i) = txPilots;
        
                % Заполнение информацией
                txMatrix(ind.data, i) = txSymbols(:, i);
            end
        
        % Обратное дискретное преобразование Френеля
        if params.style(st) == "ofdm"
            txSignal = ifft(ifftshift(txMatrix), params.N);
        else
            txSignal = params.IDFnT_matrix * txMatrix;
        end
            
        txSignalWithPrefix = zeros(params.cpLen + params.N, size(txSignal, 2));
        % Защитный интервал 
         if(prefix == "CP")
            txSignalWithPrefix(1 : params.cpLen, :) = txSignal(end - params.cpLen + 1 : end, :);
         end
         txSignalWithPrefix(params.cpLen + 1 : end, :) = txSignal;
         txSignal = txSignalWithPrefix;
    
        clear txSignalWithPrefix;
        
        % Канал связи
            % Релеевский канал (EVA)
            if(rayCh == 1)
    %             EVA.PDelaysNs = [0, 30] * 1e-9;
    %             EVA.APGdB = [0,-10];
                EVA.PDelaysNs = [0, 30, 150, 310, 370, 710, 1090, 1730, 2510] * 1e-9;
                EVA.APGdB = [0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12, -16.9];
                
                EVA.MDSHz = 0;
                rayleighChan = comm.RayleighChannel("SampleRate", params.Fs, "PathDelays", EVA.PDelaysNs, "AveragePathGains", EVA.APGdB, "MaximumDopplerShift", EVA.MDSHz);
                CR = rayleighChan([1; zeros(49, 1)]);
                rayleighChan.release();
                [~, index] = max(CR);
                rxSignal = rayleighChan([txSignal; zeros(index,1)]);
                rxSignalBeforeAWGN = rxSignal(index:index+params.N+params.cpLen-1);
                tmp_sig = [tmp_sig; rxSignalBeforeAWGN];
            else
                rxSignalBeforeAWGN = txSignal;
            end
            % АБГШ
         
            
            N0 = Eb/10^(SNR/10);
            noise = normrnd(0,sqrt(N0/2),length(rxSignalBeforeAWGN),2) * [1; 1i];
            rxSignal = rxSignalBeforeAWGN + noise;
            if params.FreqShift ~= 0
                rxSignal = rxSignal.*exp(1j*(params.FreqShift)*2*pi*(1:length(rxSignal)))';
            end
          
        % Приёмная часть
        % Обработка принятого сигнала
            % Разбиение по кадрам
            rxMatrix = reshape(rxSignal, params.N + params.cpLen, []);
            rxMatrix = rxMatrix(params.cpLen + 1 : end, :);
            % Прямое Преобразование Френеля
            if params.style(st) == "ofdm"
                rxMatrix = fftshift(fft(rxMatrix,  params.N));
            else
                rxMatrix = params.DFnT_matrix * rxMatrix * params.N;
            end
    
        
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
            txPilots1D = txPilots(:);
            H =  rxPilots1D ./ txPilots1D;
            Heq = interp1(ind.pilots, H, ind.data, 'spline');
            rxSymbols = rxSymbols ./ Heq;
    %         rxPilots = rxPilots ./ H;
        end
        % Демодуляция
            rxBitsExtended = qamdemod(rxSymbols, params.M, 'gray', 'UnitAveragePower', true, 'OutputType', 'bit');
            rxBitsExtended = rxBitsExtended(:);
        
            rxBits= rxBitsExtended(1:length(txBits));
        
        % Вычисление BER
            BER(s) = BER(s) + sum((txBits ~= rxBits));
        end
    end
    
    BER = BER/(length(txBits)*NumFrame);
    
    save(params.SaveFileName, 'h2dB', 'BER', 'params');
end