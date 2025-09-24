clear, clc

frmt = 1; % 0 – одна, 1 – много

t = linspace(-5.2,5.2, 10000);
fs = 1/length(t);
y1 = sinc(t);
y2 = sinc(t + 1.9992);
y3 = sinc(t + 2*1.9992);
y4 = sinc(t - 1.9992);
y5 = sinc(t - 2*1.9992);
figure
switch frmt
    case 1
        plot(t,y1, 'black', t, y2, 'black', t, y3, 'black', t, y4, 'black', t, y5, 'black');
    case 0
        plot(t,y1, 'black');
end
hold on
long_t = linspace(-5.5,6);
plot(long_t, zeros(length(long_t)), 'black')
yticks([]);
yticklabels([]);
xticks([]);
xticklabels([]);

% clear; clc;
% 
% N = 2;                % число поднесущих
% T = 1;                % длительность (нормируем)
% t = linspace(-5,5,1000); % временной вектор
% 
% % Формируем поднесущие sinc с ортогональными сдвигами
% figure; hold on;
% colors = lines(N);
% 
% for k = 0:N-1
%     % Смещение аргумента sinc на целое число
%     y = sinc(t - k);
%     plot(t, y, 'Color', colors(k+1,:), 'LineWidth', 1.5);
% end
% 
% % Оформление графика
% plot([-5 5], [0 0], 'k--'); % ось X
% yticks([]);
% xticks(-5:1:5);
% xlabel('Time');
% title('Orthogonal sinc subcarriers');
% 
% hold off;
% grid on;
