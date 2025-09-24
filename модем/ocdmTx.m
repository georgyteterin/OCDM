clear, clc

sdr = findsdru();

tx = comm.SDRuTransmitter(...
              'Platform',sdr.Platform, ...
              'SerialNum',sdr.SerialNum, ...
              'CenterFrequency',2.5e9, ...
              'InterpolationFactor',256);
%%
N = 1024;
cpLen = 0.25*N;
Fs = 10e6;
mapperType = 'qam'; %'psk', 'qam'
M = 4;
log2M = log2(M);
numBits = N*log2M;

mtrx = DFnTmtrx(N);
imtrx = ctranspose(mtrx);

switch mapperType
    case 'qam'
        mapper = comm.PSKModulator(M, 0, "BitInput",true);
    case 'psk'
        mapper = comm.RectangularQAMModulator(M, "BitInput",true);
end


txBits = randi([0 1], numBits, 1);
txSyms = mapper(txBits);
txSignal_ocdm = imtrx*txSyms;