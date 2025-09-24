%v1
clear, clc
semilogy(1:10, berawgn(1:10, 'psk', 4, 'nondiff'), 'blue')
hold on
semilogy([2, 3, 4, 5, 6, 7]+2.5, [0.0306, 0.0124, 0.0037, 0.0009, 1.6e-4, 2.34e-5], 'ro-')
grid on
xl = xlabel('Eb/N0, дБ');
yl = ylabel('Вероятность битовой ошибки');
xl.FontName = 'Times New Roman';
yl.FontName = 'Times New Roman';
l = legend('Теория', 'Эксперимент');
l.FontName = 'Times New Roman';
l.FontSize = 12;

%% v2
clear, clc
semilogy(1:10, berawgn(1:10, 'psk', 4, 'nondiff'), 'blue')
hold on
semilogy([2, 3, 4, 5, 6, 7, 8, 9], [0.0858, 0.069, 0.05, 0.03, 0.0124, 0.004, 0.0009, 1.15e-4], 'ro-')
grid on
xl = xlabel('Eb/N0, дБ');
yl = ylabel('Вероятность битовой ошибки');
xl.FontName = 'Times New Roman';
yl.FontName = 'Times New Roman';
l = legend('Теория', 'Эксперимент');
l.FontName = 'Times New Roman';
l.FontSize = 12;