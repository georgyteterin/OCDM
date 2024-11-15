clear, clc

N = 2^2;
input = randi([0 1], N, 1);
iniPh = exp(-1j*(pi/4));
sqrtN = sqrt(N);


% через DFT

% Teta1
Teta1 = zeros(N);
for m=1:N
    for n=1:N
        if n==m
            if mod(N, 2) == 0
                Teta1(m, n) = iniPh*exp(1j*(pi/N)*m^2);
            else
                Teta1(m, n) = iniPh*exp(1j*(pi/(4*N)))*exp(1j*(pi/N)*(m^2 + m));
            end
        end
    end
end

% Teta2
Teta2 = zeros(N);
for m=1:N
    for n=1:N
        if m==n
            if mod(N, 2) == 0
                Teta2(m, n) = exp(1j*(pi/N)*n^2);
            else
                Teta2(m, n) = exp(1j*(pi/N)*(n^2 - n));
            end
        end
    end
end

% dft matrix

dftMat = zeros(N);
for m=1:N
    for n=1:N
        dftMat(n, m) = (1/sqrtN)*exp(-1j*(2*pi/N)*m*n);
    end
end

%%
% матрицы
out1 = DFnTmtrx(input);
out2 = Teta1*dftMat*Teta2;

% преобразование

output1 = DFnTmtrx(input)*input;
tmp1 = input'*Teta1;
tmp2 = tmp1*dftMat;
output2 = tmp2*Teta2;
% output2 = Teta1*dftMat*Teta2*input;

%% без матриц

for m=1:N
    if mod(N, 2) == 0
        tmp1(m) = input(m)*iniPh*exp(1j*(pi/N)*m^2);
    else
        tmp1(m) = input(m)*iniPh*exp(1j*(pi/(4*N)))*exp(1j*(pi/N)*(m^2 + m));
    end
end

tmp2 = (1/sqrtN)*(fftshift(fft(tmp1)));

%%

f1 = (1/sqrtN)*(dftmtx(N));
f2 = dftMat;

z = 1:N;
o1 = z*f2;
o2 = sqrtN*ifft(z, N);