clear, clc

frmt_d = 16;
switch frmt_d 
    case 8
        max_size = 120;
    case 16
        max_size = 2e4;
end

frmt = "int"+frmt_d;

N = 64;
cpLen = 0.25*N;
mapperType = 'psk'; %'psk', 'qam'
M = 4;
log2M = log2(M);
numBits = N*log2M;
numFrames = 100;

mtrx = DFnTmtrx(N);
imtrx = ctranspose(mtrx);

switch mapperType
    case 'psk'
        mapper = comm.PSKModulator(M, 0, "BitInput",true);
    case 'qam'
        mapper = comm.RectangularQAMModulator(M, "BitInput",true);
end
all_data = zeros(numFrames, numBits);

% формирование преамбулы на основе М-последовательности
preambule_syms = mseq(2, 10);
preambule = upfirdn(preambule_syms, rcosdesign(0.5,16, 10), 10, 1);
% preambule = repelem(preambule_syms, 10);
% preambule = preambule/max(abs(preambule));
preambule = preambule/max(abs(preambule))*max_size;

signal = [];
pilot_bits = randi([0 1], numBits, 1);
pilot_syms = mapper(pilot_bits);
pilot_ocdm_symbol = imtrx*pilot_syms;
tmp1 = [pilot_ocdm_symbol(end - cpLen+1:end); pilot_ocdm_symbol];
pilot_ocdm_symbol = tmp1;
clear tmp1;
signal = [signal; pilot_ocdm_symbol];

for i=1:numFrames
    bits = randi([0 1], numBits, 1);
    all_data(i, :) = bits;
    syms = mapper(bits);
    ocdm_symbol = imtrx*syms;
    tmp1 = [ocdm_symbol(end - cpLen+1:end); ocdm_symbol];
    ocdm_symbol = tmp1;
    clear tmp1;
    signal = [signal; ocdm_symbol];
end

save("txData.mat", "all_data", "preambule", "pilot_ocdm_symbol", ...
    "pilot_syms", "pilot_bits")

% подготовка к записи
signal = resample(signal,2,1);
signal = signal/max(abs(signal)); % нормирование
signal = signal*max_size; % умножение на максимум
signal = [preambule; signal];

signal = repmat(signal, 10, 1);
s = zeros(2*size(signal, 1), 1); % заготовка
switch frmt_d 
    case 8
        s(1:2:end, 1) = int8(real(signal)); % синфазная составляющая
        s(2:2:end, 1) = int8(imag(signal)); % квадратурная составляющая
    case 16
        s(1:2:end) = int16(real(signal)); % синфазная составляющая
        s(2:2:end) = int16(imag(signal)); % квадратурная составляющая
end

% s = int16(s);
[fd, err]= fopen('C:\Users\гоша\Documents\VKR\tx_record.bin', 'wb');
fwrite(fd, s, frmt);
fclose(fd);

