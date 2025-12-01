classdef ChannelAnalyzer
    % Класс для расчета энергетики канала
    %
    % ОПИСАНИЕ:
    % Класс позволяет рассчитывать SNR для канала спутник-Земля 
    % без потерь в атмосфере.
    
    properties
        % Физические константы
        lightSpeedMps = 3e8;             % Скорость света [м/с]
        boltzmanConstant = 1.38e-23;     % Постоянная Больцмана [Дж/К]
        earthRadiusMeters = 6371e3;      % Радиус Земли [м]
        
        % Параметры орбиты
        orbitHeightMeters                % Высота орбиты [м]
        
        % Параметры радиолинии
        carrierFreqHz                    % Несущая частота [Гц]
        transmitPowerWatts               % Мощность передатчика [Вт]
        tGainMaxDb                       % Коэффициент усиления главного лепестка ДН [дБ]
        rTempKelvin                      % Шумовая температура приемника [К]
        receiverGainDb                   % Коэффициент усиления приемной антенны [дБ]
    end
    
    methods
        function this = ChannelAnalyzer(orbitHeightKm, carrierFreqGHz, ...
                tPowerWatts, tGainDb, rTempK, rGainDb)
            % Конструктор класса ChannelAnalyzer
            % Инициализирует и переводит параметры
            
            this.orbitHeightMeters = orbitHeightKm * 1e3;
            this.carrierFreqHz = carrierFreqGHz * 1e9;
            this.transmitPowerWatts = tPowerWatts;
            this.tGainMaxDb = tGainDb;
            this.rTempKelvin = rTempK;
            this.receiverGainDb = rGainDb; 
        end
        
        function results = calcLink(this, elevationStartDeg, elevationEndDeg, bandwidthHz)
            % calcLink - Расчет SNR rfyfkf для диапазона углов места
            %
            % ВХОДНЫЕ АРГУМЕНТЫ:
            % elevationStartDeg - Начальный угол места [град]
            % elevationEndDeg   - Конечный угол места [град]
            % bandwidthHz       - Полоса пропускания [Гц]
            %
            % ВЫХОДНЫЕ ЗНАЧЕНИЯ:
            % results - Структура с полями:
            %   .slantRangeKm - расстояние [км]
            %   .snrDb        - Отношение сигнал/шум [дБ]
            
            % Формирование вектора углов места
            pointCount = 100;
            elevationAngleDeg = linspace(elevationStartDeg, elevationEndDeg, pointCount);
            elevationAngleRad = deg2rad(elevationAngleDeg);
            
            % Расчет геометрии
            [slantRangeMeters, scanAngleRad] = this.calcGeometry(elevationAngleRad);
            
            % Расчет SNR
            snrLinear = this.calcSnrValues(slantRangeMeters, scanAngleRad, bandwidthHz);
            
            results.slantRangeKm = slantRangeMeters / 1e3;
            results.snrDb = pow2db(snrLinear);
        end
    end
    
    methods (Access = private)
        function [slantRangeMeters, nadirAngleRad] = calcGeometry(this, elevationAngleRad)
            % calcGeometry - Расчет наклонной дальности и угла сканирования
            
            orbitRadiusMeters = this.earthRadiusMeters + this.orbitHeightMeters;
            
            % Теорема синусов для расчета угла сканирования
            sinNadirAngle = (this.earthRadiusMeters .* cos(elevationAngleRad)) / orbitRadiusMeters;
            nadirAngleRad = asin(sinNadirAngle);
            
            % Центральный угол Земли
            earthCenterAngleRad = pi - (pi / 2 + elevationAngleRad) - nadirAngleRad;
            
            % Теорема косинусов для наклонной дальности
            slantRangeMeters = sqrt(this.earthRadiusMeters^2 + orbitRadiusMeters^2 ...
                -2 * this.earthRadiusMeters * orbitRadiusMeters .* cos(earthCenterAngleRad));
        end
        
        function snrLinear = calcSnrValues(this, distanceMeters, scanAngleRad, bandwidthHz)
            % calcSnrValues - Метод для расчета SNR
            
            wavelengthMeters = this.lightSpeedMps / this.carrierFreqHz;
            
            tGainMaxLinear = db2pow(this.tGainMaxDb);
            rGainLinear = db2pow(this.receiverGainDb);
            
            % Учет эффективной площади антенны планарной АФАР
            tGainEffective = tGainMaxLinear .* cos(scanAngleRad);
            
            % Потери в свободном пространстве (FSPL)
            freeSpacePathLoss = ((4 * pi * distanceMeters) ./ wavelengthMeters).^2;
            
            % Мощность теплового шума (СПМ шума на полосу)
            noisePowerWatts = this.boltzmanConstant * this.rTempKelvin * bandwidthHz;
            
            % Мощность принятого сигнала
            receivedPowerWatts = this.transmitPowerWatts .* tGainEffective .* ...
                rGainLinear ./ freeSpacePathLoss;
            
            % Итоговый SNR
            snrLinear = receivedPowerWatts ./ noisePowerWatts;
        end
    end
end