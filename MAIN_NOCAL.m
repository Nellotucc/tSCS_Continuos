clear all
close all
clc

data_dir = 'DATA\NelloL1_left'; 
mat_files = dir(fullfile(data_dir, '*.mat'));
emg_files = {mat_files.name};
emg_files = fullfile(data_dir, emg_files);

for i = 1:length(emg_files)
    fprintf('%s\n', emg_files{i});
end

%  DELSYS IMU AND EMG FREQUENCY (FIXED)
fs_IMU = 148.1481475830078;
fs_EMG = 1259.259277343750;

%% 

% Chiamata alla nuova funzione per plottare i segnali EMG per ogni muscolo
[max_values_R] = Fullgait_NoCal_subplot_R(emg_files, fs_IMU, fs_EMG);

 %IF YOU USE SUBPLOT PUT "FIGURE;" OTHERWISE DELETE IT!
[max_values_L]= Fullgait_NoCal_subplot_L(emg_files, fs_IMU, fs_EMG);


%Calculation of Gait and Stance duration
[stance_duration_R, gait_duration_R,meanstance_R,meanPercentagestance_R,meanPercentageswing_R,percentage_matrix_stanceR,percentage_matrix_swingR]=Gait_duration_R(emg_files,fs_IMU,fs_EMG);
[stance_duration_L, gait_duration_L,meanstance_L,meanPercentagestance_L,meanPercentageswing_L,percentage_matrix_stanceL,percentage_matrix_swingL]=Gait_duration_L(emg_files,fs_IMU,fs_EMG);
[ROM_hip_R,ROM_knee_R,ROM_knee_all_R,ROM_hip_all_R]= ROM_calculation_R(emg_files,fs_IMU);
[ROM_hip_L,ROM_knee_L,ROM_knee_all_L,ROM_hip_all_L]= ROM_calculation_L(emg_files,fs_IMU);



%%
%PLOT SIGNALS OF RIGHT LEG
[max_values_Rstance]= Stance_NoCal_subplot_R(emg_files, fs_IMU, fs_EMG);

% PLOT SIGNALS OF THE LEFT LEG
[max_values_Rswing]= Swing_NoCal_subplot_R(emg_files, fs_IMU, fs_EMG);

%PLOT SIGNALS OF RIGHT LEG
[max_values_Lstance]= Stance_NoCal_subplot_L(emg_files, fs_IMU, fs_EMG);

% PLOT SIGNALS OF THE LEFT LEG
[max_values_Lswing]= Swing_NoCal_subplot_L(emg_files, fs_IMU, fs_EMG);

%% FREQUENCY PLOTS
% Nomi dei muscoli
muscles_L = {'TA_L', 'VM_L', 'VL_L','RF_L', 'GM_L', 'ST_L', 'GASTRO_L', 'SO_L'};
muscles_R = {'TA_R', 'VM_R', 'VL_R', 'RF_R', 'GM_R', 'ST_R', 'GASTRO_R', 'SO_R'};

frequencies = [100, 15, 30, 50, 5, 65, 75, 0];

% Rimuovere l'undicesima colonna dai dati
max_values_L = max_values_L(:, 1:8);
max_values_R = max_values_R(:, 1:8);

% Rimuovere l'undicesima colonna dalle frequenze
frequencies = frequencies(1:8);

% Ordinare le frequenze e i corrispondenti valori di picco
[sorted_frequencies, sort_idx] = sort(frequencies);

max_values_L_sorted = max_values_L(:, sort_idx);
max_values_R_sorted = max_values_R(:, sort_idx);

% Parametri per il filtro di Savitzky-Golay
sgolay_order = 1; % Ordine del polinomio
sgolay_frame = 3; % Dimensione della finestra (più grande per un smoothing maggiore)

% Creazione delle figure per il plotting
figure;

% Plot per i muscoli della gamba sinistra
for i = 1:8
    subplot(2, 4, i); % Creazione di subplot 2x4
    
    % Applica il filtro di Savitzky-Golay
    smoothed_values_L = sgolayfilt(max_values_L_sorted(i, :), sgolay_order, sgolay_frame);
    
    semilogx(sorted_frequencies,smoothed_values_L, '-o', 'LineWidth', 2); % Plot logaritmico
    hold on;
    xline(100, '--r', 'LineWidth', 1.5); % Linea verticale tratteggiata rossa a 100 Hz
    title(muscles_L{i});
    xlabel('Frequency (Hz)');
    ylabel('Peak Value');
    grid on;
    hold off;
end

figure;

% Plot per i muscoli della gamba destra
for i = 1:8
    subplot(2, 4, i); % Creazione di subplot 2x4
    
    % Applica il filtro di Savitzky-Golay
    smoothed_values_R = sgolayfilt(max_values_R_sorted(i, :), sgolay_order, sgolay_frame);
    
    semilogx(sorted_frequencies, smoothed_values_R, '-o', 'LineWidth', 2); % Plot logaritmico
    hold on;
    xline(100, '--r', 'LineWidth', 1.5); % Linea verticale tratteggiata rossa a 100 Hz
    title(muscles_R{i});
    xlabel('Frequency (Hz)');
    ylabel('Peak Value');
    grid on;
    hold off;
end



%% STANCE/SWING RIGHT FREQUENCY PLOTS
% Nomi dei muscoli
muscles_R = {'TA_R', 'VM_R','VL_R', 'RF_R', 'GM_R', 'ST_R', 'GASTRO_R', 'SO_R'};

frequencies = [100, 15, 30, 50, 5, 65, 75, 0];

% Rimuovere l'undicesima colonna dai dati
max_values_Rstance = max_values_Rstance(:,1:8);
max_values_Rswing =  max_values_Rswing(:,1:8);

% Rimuovere l'undicesima colonna dalle frequenze
frequencies = frequencies(1:8);

% Ordinare le frequenze e i corrispondenti valori di picco
[sorted_frequencies, sort_idx] = sort(frequencies);

max_values_Rstance_sorted = max_values_Rstance(:, sort_idx);
max_values_Rswing_sorted = max_values_Rswing(:, sort_idx);

max_values_Rstance_reduced=max_values_Rstance_sorted;
max_values_Rswing_reduced=max_values_Rswing_sorted;
sorted_frequencies_reduced=sorted_frequencies;
% Parametri per il filtro di Savitzky-Golay
sgolay_order = 1; % Ordine del polinomio
sgolay_frame = 3; % Dimensione della finestra (più grande per un smoothing maggiore)

% Creazione delle figure per il plotting
figure;

% Plot per i muscoli della gamba sinistra
for i = 1:8
    subplot(2, 4, i); % Creazione di subplot 2x4
    
    % Applica il filtro di Savitzky-Golay
    smoothed_values_stance = sgolayfilt(max_values_Rstance_reduced(i, :), sgolay_order, sgolay_frame);
    
    plot(sorted_frequencies_reduced,max_values_Rstance_reduced(i, :), '-o', 'LineWidth', 2); % Plot logaritmico
    hold on;
%     xline(100, '--r', 'LineWidth', 1.5); % Linea verticale tratteggiata rossa a 100 Hz
    title('Stance phase', muscles_R{i});
    xlabel('Frequency (Hz)');
    ylabel('Peak Value');
    grid on;
    hold off;
end

figure;

% Plot per i muscoli della gamba destra
for i = 1:8
    subplot(2, 4, i); % Creazione di subplot 2x4
    
    % Applica il filtro di Savitzky-Golay
    smoothed_values_swing = sgolayfilt(max_values_Rswing_reduced(i, :), sgolay_order, sgolay_frame);
    plot(sorted_frequencies_reduced, max_values_Rswing_reduced(i, :), '-o', 'LineWidth', 2); % Plot logaritmico
    hold on;
%     xline(100, '--r', 'LineWidth', 1.5); % Linea verticale tratteggiata rossa a 100 Hz
    title('Swing phase',muscles_R{i});
    xlabel('Frequency (Hz)');
    ylabel('Peak Value');
    grid on;
    hold off;
end


%% RIGHT PERCENTAGE

frequenze_interessanti = [5, 15, 30, 50, 75, 100];
Firstvalues_stanceR=max_values_Rstance_reduced(:,1);
Firstvalues_swingR=max_values_Rswing_reduced(:,1);

Othervalues_stanceR=max_values_Rstance_reduced(:,2:end);
Othervalues_swingR=max_values_Rswing_reduced(:,2:end);


% Calcola le percentuali rispetto alla frequenza 0
percentuali_stanceR = ((Othervalues_stanceR-Firstvalues_stanceR) ./ (Firstvalues_stanceR))*100;
percentuali_swingR = ((Othervalues_swingR-Firstvalues_swingR)./ (Firstvalues_swingR))*100;

% Creazione delle figure per il plotting
figure;

% Plot per i muscoli della gamba destra (fase stance)
for i = 1:8
    subplot(2, 4, i); % Creazione di subplot 2x4
    plot(frequenze_interessanti, percentuali_stanceR(i, :), '-o', 'LineWidth', 2); % Plot percentuale
    hold on;
    title(['Stance phase - ' muscles_R{i}]);
    xlabel('Frequency (Hz)');
    ylabel('Percentage of EMG value at 0 Hz (%)');
    grid on;
    hold off;
end

figure;

% Plot per i muscoli della gamba destra (fase swing)
for i = 1:8
    subplot(2, 4, i); % Creazione di subplot 2x4
    plot(frequenze_interessanti, percentuali_swingR(i, :), '-o', 'LineWidth', 2); % Plot percentuale
    hold on;
    title(['Swing phase - ' muscles_R{i}]);
    xlabel('Frequency (Hz)');
    ylabel('Percentage of EMG value at 0 Hz (%)');
    grid on;
    hold off;
end


%% STANCE/SWING LEFT FREQUENCY PLOTS
% Nomi dei muscoli
muscles_L = {'TA_L', 'VM_L','VL_L' 'RF_L', 'GM_L', 'ST_L', 'GASTRO_L', 'SO_L'};

frequencies = [100, 15, 30, 50, 5, 65, 75, 0];

% Rimuovere l'undicesima colonna dai dati
max_values_Lstance = max_values_Lstance(:, 1:8);
max_values_Lswing = max_values_Lswing(:, 1:8);

% Rimuovere l'undicesima colonna dalle frequenze
frequencies = frequencies(1:8);

% Ordinare le frequenze e i corrispondenti valori di picco
[sorted_frequencies, sort_idx] = sort(frequencies);

max_values_Lstance_sorted = max_values_Lstance(:, sort_idx);
max_values_Lswing_sorted =  max_values_Lswing(:, sort_idx);

max_values_Lstance_reduced=max_values_Rstance_sorted;
max_values_Lswing_reduced=max_values_Rswing_sorted;
sorted_frequencies_reduced=sorted_frequencies;

% Parametri per il filtro di Savitzky-Golay
sgolay_order = 1; % Ordine del polinomio
sgolay_frame = 3; % Dimensione della finestra (più grande per un smoothing maggiore)

% Creazione delle figure per il plotting
figure;

% Plot per i muscoli della gamba sinistra
for i = 1:8
    subplot(2, 4, i); % Creazione di subplot 2x4
    
    % Applica il filtro di Savitzky-Golay
    smoothed_values_stance = sgolayfilt(max_values_Lstance_sorted(i, :), sgolay_order, sgolay_frame);
    
    plot(sorted_frequencies_reduced,max_values_Lstance_reduced(i, :), '-o', 'LineWidth', 2); % Plot logaritmico
    hold on;
%     xline(100, '--r', 'LineWidth', 1.5); % Linea verticale tratteggiata rossa a 100 Hz
    title('Stance phase', muscles_L{i});
    xlabel('Frequency (Hz)');
    ylabel('Peak Value');
    grid on;
    hold off;
end

figure;

% Plot per i muscoli della gamba destra
for i = 1:8
    subplot(2, 4, i); % Creazione di subplot 2x4
    
    % Applica il filtro di Savitzky-Golay
    smoothed_values_swing = sgolayfilt(max_values_Lswing_sorted(i, :), sgolay_order, sgolay_frame);
    
    plot(sorted_frequencies_reduced, max_values_Lswing_reduced(i, :), '-o', 'LineWidth', 2); % Plot logaritmico
    hold on;
%     xline(100, '--r', 'LineWidth', 1.5); % Linea verticale tratteggiata rossa a 100 Hz
    title('Swing phase',muscles_L{i});
    xlabel('Frequency (Hz)');
    ylabel('Peak Value');
    grid on;
    hold off;
end

%% LEFT IN PERCENTAGE
% Nomi dei muscoli
muscles_L = {'TA_L', 'VM_L','VL_L', 'RF_L', 'GM_L', 'ST_L', 'GASTRO_L', 'SO_L'};

% Frequenze da considerare per il plot come percentuali
frequenze_interessanti = [5, 15, 30, 50, 75, 100];

Firstvalues_stanceL=max_values_Lstance_reduced(:,1);
Firstvalues_swingL=max_values_Lstance_reduced(:,1);

Othervalues_stanceL=max_values_Lstance_reduced(:,2:end);
Othervalues_swingL=max_values_Lswing_reduced(:,2:end);


% Calcola le percentuali rispetto alla frequenza 0
percentuali_stanceL = ((Othervalues_stanceL-Firstvalues_stanceL) ./ (Firstvalues_stanceL))*100;
percentuali_swingL = ((Othervalues_swingL-Firstvalues_swingL)./ (Firstvalues_swingL))*100;


% Creazione delle figure per il plotting
figure;

% Plot per i muscoli della gamba destra (fase stance)
for i = 1:8
    subplot(2, 4, i); % Creazione di subplot 2x4
    plot(frequenze_interessanti, percentuali_stanceL(i, :), '-o', 'LineWidth', 2); % Plot percentuale
    hold on;
    title(['Stance phase - ' muscles_L{i}]);
    xlabel('Frequency (Hz)');
    ylabel('Percentage of EMG value at 0 Hz (%)');
    grid on;
    hold off;
end

figure;

% Plot per i muscoli della gamba destra (fase swing)
for i = 1:8
    subplot(2, 4, i); % Creazione di subplot 2x4
    plot(frequenze_interessanti, percentuali_swingL(i, :), '-o', 'LineWidth', 2); % Plot percentuale
    hold on;
    title(['Swing phase - ' muscles_L{i}]);
    xlabel('Frequency (Hz)');
    ylabel('Percentage of EMG value at 0 Hz (%)');
    grid on;
    hold off;
end
