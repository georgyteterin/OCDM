figure();
semilogy(2:9, berawgn(2:9, 'psk', 4, 'nondiff'), 'blue')
hold on
semilogy([2 3 4 4.5 5 5.5 6]+2.5, [0.0343 0.0146 0.0033 0.0016 7.8e-4 3.029e-4 9.4e-5], 'ro-')
grid on
xl = xlabel('Eb/N0, дБ');
xl.FontName = "Times New Roman";
yl = ylabel("Вероятность битовой ошибки");
l = legend('Теория', 'Эксперимент');
l.FontName = "Times New Roman";
l.FontSize = 12;
