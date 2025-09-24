clear, clc

% Параметры
N = 16;          % количество поднесущих
T = 1;           % период сигнала
fs = 1000;       % частота дискретизации
t = linspace(0, T, fs); % временная ось

% --- Генерация OFDM поднесущих (синусоиды с разной частотой) ---
ofdm_subcarriers = zeros(N, length(t));
for k = 1:N
    freq = k/ T; % частота k-й поднесущей
    ofdm_subcarriers(k, :) = sin(2*pi*freq*t);
end

% --- Генерация OCDM поднесущих ---
% ЛЧМ сигнал с разной начальной фазой для каждой поднесущей
ocdm_subcarriers = zeros(N, length(t));
B = N / T; % ширина полосы ЛЧМ
for k = 1:N+1
    phi0 = 16*pi*(k-1)/N; % начальная фаза для k-й поднесущей
    % ЛЧМ сигнал: s(t) = exp(j*pi*B*t^2/T + j*phi0)
    % Для визуализации берем реальную часть
%     ocdm_subcarriers(k, :) = cos(pi*B*t.^2/T + phi0);
    ocdm_subcarriers(k, :) = exp(1j*pi/4)*exp(-1j*pi*(N/T^2)*(t-(k-1)*(T/N)).^2);
end

%% Отдельно OCDM

figure();

hold on;
for k = 2:N+1
    y_ind(k) = 3*(k-1);
    plot(t, real(ocdm_subcarriers(k,:)) - y_ind(k), 'b'); % сдвигаем по оси y для наглядности
end
xl = xlabel('Время');
xl.FontName = 'Times New Roman';
xl.FontSize = 12;
yl = ylabel('ЛЧМ-сигналы');
yl.FontName = 'Times New Roman';
yl.FontSize = 12;
yticks([]);
yticklabels([]);
% ylim([-1 N-5]);
xticks([0, 0.5, 1]);
xticklabels({'0', 'T/2', 'T'});

hold off;


%% --- Построение графиков ---

% figure;

% % Левая часть: OFDM поднесущие
subplot(1,2,1);
hold on;
index = [1 2 3 16];
for k = 1:length(index)
    y_ind(k) = 3*(k-1);
    if index(k) ~=  index(end-1)
        plot(t, ofdm_subcarriers(index(k),:) + 3*(k-1), 'b'); % сдвигаем по оси y для наглядности
    else
        plot(t, NaN)
    end
end
xl = xlabel('Время');
xl.FontName = 'Times New Roman';
xl.FontSize = 12;
yl = ylabel('Поднесущая');
yl.FontName = 'Times New Roman';
yl.FontSize = 12;
yticks(y_ind);
yticklabels({'#1', '#2', '...', '#16'});
ylim([-1 N-5]);
xticks([0, 0.5, 1]);
xticklabels({'0', 'T/2', 'T'});

hold off;

% % Правая часть: OCDM поднесущие
subplot(1,2,2);
hold on;
index = [1 2 3 17];
for k = 1:length(index)
    y_ind(k) = 3*(k-1);
    if index(k) ~= index(end-1)
        plot(t, real(ocdm_subcarriers(index(k),:)) + 3*(k-1), 'b'); % сдвигаем по оси y для наглядности
    else
        plot(t, NaN)
    end
end
xl = xlabel('Время');
xl.FontName = 'Times New Roman';
xl.FontSize = 12;
yl = ylabel('Поднесущая');
yl.FontName = 'Times New Roman';
yl.FontSize = 12;
yticks(y_ind);
yticklabels({'#1', '#2', '...', '#16'});
ylim([-1 N-5]);
xticks([0, 0.5, 1]);
xticklabels({'0', 'T/2', 'T'});
hold off;
