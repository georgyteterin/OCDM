function DrawBER(DirName, style)

    arguments
        DirName = 'DefaultResults';
        style = "clear";
    end

    % Получим информацию о содержимом директории
    Listing = dir(DirName);

    % Инициализируем cell-массив для хранения имён файлов, из которых потом
    % сделаем легенду
    Names = cell(0);

    % Цикл по всем файлам директории
    for k = 1:length(Listing)
    % Надо проверять, чтобы рассматриваемый элемент был файлом и
    % имел расширение mat
        if ~Listing(k).isdir
            FName = Listing(k).name;
            if length(FName) > 4
                if isequal(FName(end-3:end), '.mat')
                    % Добавим имя файла к списку
                        Names{end+1} = FName(1:end-4); %#ok<AGROW>
                end
            end
        end
    end

    if isempty(Names)
        error('Не найдены файлы с результатами!');
    end
    f  = cell(1, 1);
    ax = cell(1, 1);
    f{1} = figure;
    ax{1} = axes;


    % Цикл по всем уже известным файлам
    for k = 1:length(Names)
        % Загрузка результатов
        data = load([DirName, '\', Names{k}, '.mat']);
        % Прорисовка без затирания старых рисунков    
        figure(f{1});
            hold on;
            plot(data.h2dB, data.BER, ...
                'LineWidth', 1, 'MarkerSize', 8, ...
                'Marker', '.');
%             plot(data.BER, data.h2dB, ...
%                 'LineWidth', 1, 'MarkerSize', 8, ...
%                 'Marker', '.');
    end

    % Добавим сетку
    grid on;

    % Сделаем традиционный масштаб по оси ординат
        set(ax{1}, 'YScale', 'log');

        % Подпишем рисунок и ось абсцисс
        title('BER');
        xlabel('{\ith}^2 (dB)');
        
        AddNames = cell(0);
        h2dB = 0:0.1:14.4;
        BER = berawgn(h2dB, 'qam', 16);
        plot(h2dB, BER);
        AddNames{end+1} = '16-QAM'; %#ok<AGROW>

        % Добавим легенду
        legend([Names, AddNames], 'Interpreter', 'none');

%     data = load(filename);
%     BERth = berawgn(data.h2dB, "qam", data.params.M);
%     
% %     figure();
%     title('');
%     if style == "fit"
%         berfit(data.h2dB, data.BER);
%     elseif style == "clear"
%         semilogy(data.h2dB, data.BER);
%     end
%     grid on;
%     hold on;
%     semilogy(data.h2dB, BERth, '--');
%     ylim([10^-6 1]);
%     % legend AutoUpdate on;
%     if style == "fit"
%         if data.params.style == 'ocdm'
%             legend('', [ 'OCDM 16-QAM частотный сдвиг = ' num2str(data.params.FreqShift) 'Pi'], '16-QAM');
%         else
%             legend('', [ 'OFDM 16-QAM частотный сдвиг = ' num2str(data.params.FreqShift) 'Pi'], '16-QAM');
%         end
%     else
%         if data.params.style == "ocdm"
%             legend([ 'OCDM 16-QAM частотный сдвиг = ' num2str(data.params.FreqShift) 'Pi'], '16-QAM');
%         else
%             legend([ 'OFDM 16-QAM частотный сдвиг = ' num2str(data.params.FreqShift) 'Pi'], '16-QAM');
%         end
%     end
end