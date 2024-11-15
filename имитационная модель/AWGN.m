%% инициализация

clear, clc

ModulationOrder = 16;
k = log2(ModulationOrder);
EbN0 = 0:0.5:10;
FnIsTransparant = false;
ChannelIsTransparant = false;
SP = 1;
N = 20000;
subN = N/20;
subLen = N/subN;
% sub = zeros(subLen, subN);


% source

inputSyms = randi([0 ModulationOrder-1], 1, N);
inputBits = int2bit(inputSyms, log2(ModulationOrder));

%SP

if ~SP
    SourceOut = inputSyms;
else
    SourceOut = reshape(inputSyms, [], subN);
end

% mapper
if ~SP
    MapperOut = qammod(SourceOut, ModulationOrder);
else
    for i=1:subN
        MapperOut(:, i) = qammod(SourceOut(:, i), ModulationOrder);
    end
end

% IDFnt

if FnIsTransparant
    InvFresnelOut =  reshape(MapperOut, [], 1);
else
    if ~SP
        I = real(MapperOut);
        Q = imag(MapperOut);
        TransformedI = IDFnT(I); 
        TransformedQ = IDFnT(Q);
        InvFresnelOut(:, i) = TransformedI + 1j*TransformedQ;
    else
        for i=1:subN
            I = real(MapperOut(:, i));
            Q = imag(MapperOut(:, i));
            TransformedI = IDFnT(I); 
            TransformedQ = IDFnT(Q);
            InvFresnelOut(:, i) = TransformedI + 1j*TransformedQ;
%             InvFresnelOut2(:, i) = IDFnT(MapperOut(:, i));
        end
    end
    %PS
    InvFresnelOut = reshape(InvFresnelOut, 1, []);
end

% GetSomeStats 
for SNR = 1:length(EbN0)
    % Channel
    if ChannelIsTransparant 
        ChannelOut = InvFresnelOut;
    else
        constellation = unique(MapperOut);
        Esym = mean(abs(constellation).^2);
        Eb = Esym/log2(sum(ModulationOrder));
        N0 = Eb/10^(EbN0(SNR)/10);

        noise = (randn(1,length(InvFresnelOut))+1j*randn(1,length(InvFresnelOut))).*sqrt(N0/2);
        sigRx = InvFresnelOut+noise;
    end

    % SP
    if SP
        ChannelOut = reshape(sigRx, [], subN);
    end

    % DFnT
    
    if FnIsTransparant
        FresnelOut = zeros(size(sigRx));
        FresnelOut = sigRx;
    else
        FresnelOut = zeros(size(ChannelOut));
        for i=1:subN
            I = real(ChannelOut(:, i));
            Q = imag(ChannelOut(:, i));
            TransformedI = DFnT(I); 
            TransformedQ = DFnT(Q);
            FresnelOut(:, i) = TransformedI + 1j*TransformedQ;
        end
    end
    %PS
    FresnelOut = reshape(FresnelOut, 1, []);
    % demapper
    outputBits = qamdemod(FresnelOut, ModulationOrder, "OutputType", "bit");


    % statistics
%     Out = int2bit(Out, k);
    inputBits = reshape(inputBits, [], 1);
    outputBits = reshape(outputBits, [], 1);
    BER(SNR) = sum(inputBits ~= outputBits)/length(outputBits);
%     BER(SNR) = biterr(Out, Bits);
%     clear SourceOut; 
%     clear MapperOut;
%     clear InvFresnelOut;
%     clear ChannelOut;
%     clear FresnelOut;
%     clear DemapperOut; 
end


%% drawBER;
figure();
title('');
berfit(EbN0, BER)
hold on
ber_theory = berawgn(EbN0,'qam', 16);
semilogy(EbN0,ber_theory);
legend('', 'OCDM approx', 'QAM-4');

%% Spectrum


