% matlab скрипт для расчета энергетики канала
% 
% runLinkAnalysis - Скрипт запуска анализа энергетики радиолинии
%
% Описание
% Создает обьект класса расчёта канала, запускает расчет SNR для различной
% ширины полосы сигнала. Строит график SNR от угла места.

clc;
clear;
close all;
% Не успел: вынести на вход класса параметры антенны и количество точек
% моделирования, написать тест под новый модуль, нормировку на КУ заданный
% в ТЗ(Если мы считаем его истинным)
%% Инициализация параметров
orbitHeightKm = 800;    % Высота орбиты [км]
freqGHz = 12.5;         % Несущая частота [ГГц]
tPowerWatts = 50;       % Мощность передатчика [Вт]
tGainDb = 30.1;         % Коэффициент усиления главного лепестка ДН [дБ]
% Можно еще добавить проверку т.к. есть доп. параметры (антенны) pow2db(32*32)=30.1 
% И lambda=c/f >=2d
rTempKelvin = 530;      % Шумовая температура приемника [К]
rGainDb = 0;            % Коэффициент усиления приемной антенны [дБ]

% Создание обьекта класса
channelAnalyzer = ChannelAnalyzer(orbitHeightKm, freqGHz, ...
    tPowerWatts, tGainDb, rTempKelvin, rGainDb);

% Параметры моделирования
elevationStartDeg = 90;     % Начальный угол места
elevationEndDeg = 27.5;   % Конечный угол места
bandwidthListHz = [1.44, 10, 20, 50] * 1e6;  % Ширина полосы
bandwidthCount = length(bandwidthListHz);

%% Расчет и визуализация

figure('Name', 'Расчет SNR канала спутник-Земля', 'Color', 'w');
hold on;
grid on;


for bwIdx = 1:bandwidthCount
    currentBandwidthHz = bandwidthListHz(bwIdx);
    
    % Вызов метода расчета SNR
    results = channelAnalyzer.calcLink(elevationStartDeg, ...
        elevationEndDeg, currentBandwidthHz);
    
    plot(results.slantRangeKm, results.snrDb, 'LineWidth', 2, ...
        'DisplayName', sprintf('Bandwidth: %.2f MHz', currentBandwidthHz / 1e6));
end

%% Оформление графика
legend('Location', 'best');
xlim([min(results.slantRangeKm), max(results.slantRangeKm)]);
title('Зависимость SNR от наклонной дальности');
xlabel('Наклонная дальность [км]');
ylabel('Отношение сигнал/шум  [дБ]');
legend();
hold off;