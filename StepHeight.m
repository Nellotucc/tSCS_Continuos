
clear all
close all
clc


load DATA\23_05_DOME\DOME_Withoutstim2.mat;

%% Cambio delle prime 8 righe 
% Supponiamo che la tua matrice si chiami Data
% Estrai le dimensioni della matrice
[nRows, nCols] = size(Data);

% Controlla se ci sono almeno 8 righe
if nRows >= 8
    % Estrai le prime 8 righe
    firstEightRows = Data(1:8, :);
    
    % Estrai le restanti righe
    remainingRows = Data(9:end,:);
    
    % Concatenale per ottenere la nuova matrice
    newData = [remainingRows; firstEightRows];
else
    error('matrix doesnt have 8 rows');
end

% Assegna la nuova matrice a Data
Data = newData;
fs_IMU=148.1481475830078;

[contact_times, heelstrikes_times] = findHeelContactTimes_L(Data, fs_IMU);

% Numero di cicli di camminata
num_cycles = length(contact_times) - 1;

% Prealloca le variabili per velocità e posizione
velocity_segments = cell(num_cycles, 1);
position_segments = cell(num_cycles, 1);

% Itera attraverso i cicli di camminata
for i = 1:num_cycles
    t1 = contact_times(i);
    t2 = contact_times(i + 1);

    % Calcola gli indici corrispondenti
    index1 = round(t1 * fs_IMU);
    index2 = round(t2 * fs_IMU);

    l=length(Data(67,:));
    f=Data(67,1:(l/8.5));
%     pizza=f-mean(f);

    % Estrai il segmento di segnale tra gli indici calcolati
    segment = f(index1:index2);
    segment = segment*9.81;
    % Calcola il vettore del tempo corrispondente al segmento
    t_segment = (index1:index2) / fs_IMU;

    % Filtraggio del segnale (opzionale, da sostituire con i tuoi parametri di filtro)
    fc = 5; % Frequenza di taglio del filtro passa-basso (Hz)
    [b, a] = butter(2, fc / (fs_IMU / 2)); % Filtro passa-basso Butterworth di ordine 2
    filtered_segment = filtfilt(b, a, segment);

    % Controlla la direzione dell'accelerazione
    % Se necessario, inverti il segno del segnale
%     filtered_segment = -filtered_segment;

    % Integra il segnale di accelerazione per ottenere la velocità
    velocity_segment = cumtrapz(t_segment, filtered_segment);

    % Correggi la velocità per assumere che la velocità iniziale e finale sia zero
    velocity_segment = velocity_segment - mean(velocity_segment([1 end]));

    % Integra il segnale di velocità per ottenere la posizione
    position_segment = cumtrapz(t_segment, velocity_segment);

    % Memorizza i risultati
    velocity_segments{i} = velocity_segment;
    position_segments{i} = position_segment;
end

% Normalizza le posizioni per il ciclo di camminata e calcola la media
normalized_positions = zeros(num_cycles, length(linspace(0, 1, length(position_segments{1}))));
for i = 1:num_cycles
    % Interpola ogni segmento di posizione alla stessa lunghezza
    normalized_positions(i, :) = interp1(linspace(0, 1, length(position_segments{i})), position_segments{i}, linspace(0, 1, length(position_segments{1})));
end
mean_position = mean(normalized_positions);

% Plot del segnale, velocità e posizione per ogni ciclo di camminata
figure;
subplot(2, 1, 1);
hold on;
for i = 1:num_cycles
    % Normalizza il tempo per il ciclo di camminata
    cycle_time = linspace(0, 1, length(velocity_segments{i}));
    plot(cycle_time, velocity_segments{i}, 'LineWidth', 2);
end
xlabel('Ciclo di camminata');
ylabel('Velocità');
title('Velocità Integrata per ogni ciclo di camminata');
hold off;

subplot(2, 1, 2);
hold on;
for i = 1:num_cycles
    % Normalizza il tempo per il ciclo di camminata
    cycle_time = linspace(0, 1, length(position_segments{i}));
    plot(cycle_time, position_segments{i}, 'LineWidth', 2);
end
% Aggiungi la media in nero
plot(linspace(0, 1, length(mean_position)), mean_position, 'k', 'LineWidth', 3);
xlabel('Ciclo di camminata');
ylabel('Posizione');
title('Posizione Integrata per ogni ciclo di camminata');
legend([arrayfun(@(x) sprintf('Ciclo %d', x), 1:num_cycles, 'UniformOutput', false), 'Media']);
hold off;

