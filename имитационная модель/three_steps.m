clear, clc

N = 2^10;
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
        dftMat(n, m) = (1/sqrtN)*exp(-1j*(2*pi/N)*(m - 1)*(n - 1));
    end
end

DFnTmat = DFnTmtrx(input);
%%
% матрицы

out1 = input'*DFnTmat;

out2 = input'*Teta1*dftMat*Teta2;

% преобразование

% tic;output1 = DFnT(input);toc;
tic;output1 = DFnTmat*input;toc;

tic;
tmp11 = Teta1*input;
tmp12 = dftMat*tmp11;
output2 = Teta2*tmp12;toc;

tic;
tmp21 = Teta1*input;
tmp22 = fft(tmp21, N);
tmp22 = (1/sqrtN)*tmp22;
output3 = Teta2*tmp22;toc;
