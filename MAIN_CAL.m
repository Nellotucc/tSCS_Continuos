close all
clear all
clc

%DIRECTORY SPECIFICATION
data_dir = 'DATA\SUB1'; 
mat_files = dir(fullfile(data_dir, '*.mat'));
emg_files = {mat_files.name};
emg_files = fullfile(data_dir, emg_files);


for i = 1:length(emg_files)
    fprintf('%s\n', emg_files{i});
end

%  DELSYS IMU AND EMG FREQUENCY (FIXED)
fs_IMU = 148.1481475830078;
fs_EMG = 1259.259277343750; 

%% SUBPLOT + STANCE/SWING DURATION + ROM KNEE/HIP


% %PLOT SIGNALS OF RIGHT LEG
% [max_values_R]= Fullgait_Cal_subplot_R(emg_files, fs_IMU, fs_EMG);
% 
% % PLOT SIGNALS OF THE LEFT LEG
% [max_values_L]= Fullgait_Cal_subplot_L(emg_files, fs_IMU, fs_EMG);

%PLOT SIGNALS OF RIGHT LEG
[max_values_Rstance]= Stance_Cal_subplot_R(emg_files, fs_IMU, fs_EMG);

% PLOT SIGNALS OF THE LEFT LEG
[max_values_Rswing]= Swing_Cal_subplot_R(emg_files, fs_IMU, fs_EMG);

%PLOT SIGNALS OF RIGHT LEG
[max_values_Lstance]= Stance_Cal_subplot_L(emg_files, fs_IMU, fs_EMG);

% PLOT SIGNALS OF THE LEFT LEG
[max_values_Lswing]= Swing_Cal_subplot_L(emg_files, fs_IMU, fs_EMG);


%Gait(stance & swing)duration, knee & swing ROM (in stance & swing) 
[stance_duration_R,swing_duration_R, gait_duration_R,meanstance_R,meanPercentagestance_R,meanPercentageswing_R,percentage_matrix_stanceR,percentage_matrix_swingR]=Gait_duration_R(emg_files,fs_IMU,fs_EMG);
[stance_duration_L,swing_duration_L, gait_duration_L,meanstance_L,meanPercentagestance_L,meanPercentageswing_L,percentage_matrix_stanceL,percentage_matrix_swingL]=Gait_duration_L(emg_files,fs_IMU,fs_EMG);
[ROM_hip_R,ROM_knee_R,ROM_knee_all_R,ROM_hip_all_R]= ROM_calculation_R(emg_files,fs_IMU);
[ROM_hip_L,ROM_knee_L,ROM_knee_all_L,ROM_hip_all_L]= ROM_calculation_L(emg_files,fs_IMU);
[ ROM_hip_stance_L, ROM_hip_swing_L, ROM_knee_stance_L, ROM_knee_swing_L, ROM_hip_stance_all_L, ROM_hip_swing_all_L, ROM_knee_stance_all_L, ROM_knee_swing_all_L] = ROM_L(emg_files, fs_IMU);
[ ROM_hip_stance_R, ROM_hip_swing_R, ROM_knee_stance_R, ROM_knee_swing_R, ROM_hip_stance_all_R, ROM_hip_swing_all_R, ROM_knee_stance_all_R, ROM_knee_swing_all_R] = ROM_R(emg_files, fs_IMU);



meanValues_R = zeros(15,1);
meanValues_L=  zeros(15,1);

for i = 1:15
    non_zero_elements = gait_duration_R(i, gait_duration_R(i,:) ~= 0);
    if ~isempty(non_zero_elements)
        meanValues_R(i) = mean(non_zero_elements);
    else
        meanValues_R(i) = NaN; 
    end
end

for i = 1:15
    non_zero_elements = gait_duration_L(i, gait_duration_L(i,:) ~= 0);
    if ~isempty(non_zero_elements)
        meanValues_L(i) = mean(non_zero_elements);
    else
        meanValues_L(i) = NaN; 
    end
end


meanPercentageswing_L_DIST=meanPercentageswing_L([1:7,15]);
meanPercentagestance_L_DIST=meanPercentagestance_L([1:7,15]);

meanPercentageswing_L_PROX=meanPercentageswing_L(8:15,:);
meanPercentagestance_L_PROX=meanPercentagestance_L(8:15,:);


meanPercentageswing_R_DIST=meanPercentageswing_R([1:7,15]);
meanPercentagestance_R_DIST=meanPercentagestance_R([1:7,15]);
meanPercentageswing_R_PROX=meanPercentageswing_R(8:15,:);
meanPercentagestance_R_PROX=meanPercentagestance_R(8:15,:);

frequencies = [100,15,30,50,5,65,75,0];

% Sorting frequencies and relative percentage
[sorted_frequencies, sort_idx] = sort(frequencies);
meanPercentagestance_R_DISTsorted = meanPercentagestance_R_DIST(sort_idx);
meanPercentageswing_R_DISTsorted  = meanPercentageswing_R_DIST(sort_idx);

meanPercentagestance_R_PROXsorted = meanPercentagestance_R_PROX(sort_idx);
meanPercentageswing_R_PROXsorted  = meanPercentageswing_R_PROX(sort_idx);


meanPercentagestance_L_DISTsorted = meanPercentagestance_L_DIST(sort_idx);
meanPercentageswing_L_DISTsorted  = meanPercentageswing_L_DIST(sort_idx);

meanPercentagestance_L_PROXsorted = meanPercentagestance_L_PROX(sort_idx);
meanPercentageswing_L_PROXsorted  = meanPercentageswing_L_PROX(sort_idx);

figure;
x_positions = 1:length(sorted_frequencies);
bar_width = 0.3;
bar(x_positions - bar_width/2, meanPercentagestance_R_DISTsorted, bar_width, 'b');
hold on;
x_positions = 1:length(sorted_frequencies);
bar(x_positions + bar_width/2, meanPercentageswing_R_DISTsorted, bar_width, 'r');
xlabel('Frequencies');
ylabel('Percentage');
xticks(x_positions);
xticklabels(cellfun(@(x) [num2str(x) ' Hz'], num2cell(sorted_frequencies), 'UniformOutput', false));
title('Phases percentage R- DISTAL');
legend('Stance Phase', 'Swing Phase', 'Location', 'northwest'); 


figure;
x_positions = 1:length(sorted_frequencies);
bar_width = 0.3;
bar(x_positions - bar_width/2, meanPercentagestance_R_PROXsorted, bar_width, 'b');
hold on;
x_positions = 1:length(sorted_frequencies);
bar(x_positions + bar_width/2, meanPercentageswing_R_PROXsorted, bar_width, 'r');
xlabel('Frequencies');
ylabel('Percentage');
xticks(x_positions);
xticklabels(cellfun(@(x) [num2str(x) ' Hz'], num2cell(sorted_frequencies), 'UniformOutput', false));
title('Phases percentage R- PROX');
legend('Stance Phase', 'Swing Phase', 'Location', 'northwest');

figure;
x_positions = 1:length(sorted_frequencies);
bar_width = 0.3;
bar(x_positions - bar_width/2, meanPercentagestance_L_PROXsorted, bar_width, 'b');
hold on;
x_positions = 1:length(sorted_frequencies);
bar(x_positions + bar_width/2, meanPercentageswing_L_PROXsorted, bar_width, 'r');
xlabel('Frequencies');
ylabel('Percentage');
xticks(x_positions);
xticklabels(cellfun(@(x) [num2str(x) ' Hz'], num2cell(sorted_frequencies), 'UniformOutput', false));
title('Phases percentage L- PROX');
legend('Stance Phase', 'Swing Phase', 'Location', 'northwest'); 


figure;
x_positions = 1:length(sorted_frequencies);
bar_width = 0.3;
bar(x_positions - bar_width/2, meanPercentagestance_L_DISTsorted, bar_width, 'b');
hold on;
x_positions = 1:length(sorted_frequencies);
bar(x_positions + bar_width/2, meanPercentageswing_L_DISTsorted, bar_width, 'r');
xlabel('Frequencies');
ylabel('Percentage');
xticks(x_positions);
xticklabels(cellfun(@(x) [num2str(x) ' Hz'], num2cell(sorted_frequencies), 'UniformOutput', false));
title('Phases percentage L- DIST');
legend('Stance Phase', 'Swing Phase', 'Location', 'northwest'); 



ROM_knee_TOT=(ROM_knee_R+ROM_knee_L)/2;
ROM_hip_TOT=(ROM_hip_R+ROM_hip_L)/2;

ROM_knee_stance_TOT=(ROM_knee_stance_R+ROM_knee_stance_L)/2;
ROM_hip_stance_TOT=(ROM_hip_stance_R+ROM_hip_stance_L)/2;

ROM_knee_swing_TOT=(ROM_knee_swing_R+ROM_knee_swing_L)/2;
ROM_hip_swing_TOT=(ROM_hip_swing_R+ROM_hip_swing_L)/2;



ROM_knee_stance_TOTPROX=ROM_knee_stance_TOT(8:15,:);
ROM_hip_stance_TOTPROX=ROM_hip_stance_TOT(8:15,:);

ROM_knee_stance_TOTDIST=ROM_knee_stance_TOT([1:7,15]);
ROM_hip_stance_TOTDIST=ROM_hip_stance_TOT([1:7,15]);



ROM_knee_swing_TOTPROX=ROM_knee_swing_TOT(8:15,:);
ROM_hip_swing_TOTPROX=ROM_hip_swing_TOT(8:15,:);

ROM_knee_swing_TOTDIST=ROM_knee_swing_TOT([1:7,15]);
ROM_hip_swing_TOTDIST=ROM_hip_swing_TOT([1:7,15]);


ROM_knee_stance_TOTPROX_sorted=ROM_knee_stance_TOTPROX(sort_idx)';
ROM_hip_stance_TOTPROX_sorted=ROM_hip_stance_TOTPROX(sort_idx)';
ROM_knee_stance_TOTDIST_sorted=ROM_knee_stance_TOTDIST(sort_idx)';
ROM_hip_stance_TOTDIST_sorted=ROM_hip_stance_TOTDIST(sort_idx)';


ROM_knee_swing_TOTPROX_sorted=ROM_knee_swing_TOTPROX(sort_idx)';
ROM_hip_swing_TOTPROX_sorted=ROM_hip_swing_TOTPROX(sort_idx)';
ROM_knee_swing_TOTDIST_sorted=ROM_knee_swing_TOTDIST(sort_idx)';
ROM_hip_swing_TOTDIST_sorted=ROM_hip_swing_TOTDIST(sort_idx)';





 ROM_knee_R_PROX= ROM_knee_R(8:15,:);
 ROM_knee_L_PROX= ROM_knee_L(8:15,:);

 ROM_knee_R_DIST= ROM_knee_R([1:7,15]);
 ROM_knee_L_DIST= ROM_knee_L([1:7,15]);


 ROM_hip_R_PROX= ROM_hip_R(8:15,:);
 ROM_hip_L_PROX= ROM_hip_L(8:15,:);

 ROM_hip_R_DIST= ROM_hip_R([1:7,15]);
 ROM_hip_L_DIST= ROM_hip_L([1:7,15]);


 ROM_knee_PROX= ROM_knee_TOT(8:15,:);
 ROM_knee_DIST= ROM_knee_TOT([1:7,15]);
 ROM_hip_PROX= ROM_hip_TOT(8:15,:);
 ROM_hip_DIST= ROM_hip_TOT([1:7,15]);




 ROM_knee_swing_R_PROX=ROM_knee_swing_R(8:15,:);
 ROM_knee_stance_R_PROX=ROM_knee_stance_R(8:15,:);

 ROM_knee_swing_R_DIST=ROM_knee_swing_R([1:7,15]);
 ROM_knee_stance_R_DIST=ROM_knee_stance_R([1:7,15]);

 ROM_knee_swing_L_PROX=ROM_knee_swing_L(8:15,:);
 ROM_knee_stance_L_PROX=ROM_knee_stance_L(8:15,:);

 ROM_knee_swing_L_DIST=ROM_knee_swing_L([1:7,15]);
 ROM_knee_stance_L_DIST=ROM_knee_stance_L([1:7,15]); 
 

ROM_hip_stance_R_PROX=ROM_hip_stance_R(8:15,:);
ROM_hip_swing_R_PROX=ROM_hip_swing_R(8:15,:);

ROM_hip_stance_R_DIST=ROM_hip_stance_R([1:7,15]);
ROM_hip_swing_R_DIST=ROM_hip_swing_R([1:7,15]);


ROM_hip_stance_L_PROX=ROM_hip_stance_L(8:15,:);
ROM_hip_swing_L_PROX=ROM_hip_swing_L(8:15,:);

ROM_hip_stance_L_DIST=ROM_hip_stance_L([1:7,15]);
ROM_hip_swing_L_DIST=ROM_hip_swing_L([1:7,15]);



frequencies = [100,15,30,50,5,65,75,0];
[sorted_frequencies, sort_idx] = sort(frequencies);
ROM_knee_R_sortedPROX = ROM_knee_R_PROX(sort_idx);
ROM_knee_L_sortedPROX = ROM_knee_L_PROX(sort_idx);
ROM_knee_R_sortedDIST = ROM_knee_R_DIST(sort_idx);
ROM_knee_L_sortedDIST = ROM_knee_L_DIST(sort_idx);

ROM_hip_R_sortedPROX = ROM_hip_R_PROX(sort_idx);
ROM_hip_L_sortedPROX = ROM_hip_L_PROX(sort_idx);
ROM_hip_R_sortedDIST = ROM_hip_R_DIST(sort_idx);
ROM_hip_L_sortedDIST = ROM_hip_L_DIST(sort_idx);




ROM_knee_stance_RsortedPROX=ROM_knee_stance_R_PROX(sort_idx);
ROM_knee_swing_RsortedPROX=ROM_knee_swing_R_PROX(sort_idx);
ROM_knee_stance_LsortedPROX=ROM_knee_stance_L_PROX(sort_idx);
ROM_knee_swing_LsortedPROX=ROM_knee_swing_L_PROX(sort_idx);


ROM_knee_stance_RsortedDIST=ROM_knee_stance_R_DIST(sort_idx);
ROM_knee_swing_RsortedDIST=ROM_knee_swing_R_DIST(sort_idx);
ROM_knee_stance_LsortedDIST=ROM_knee_stance_L_DIST(sort_idx);
ROM_knee_swing_LsortedDIST=ROM_knee_swing_L_DIST(sort_idx);



ROM_hip_stance_RsortedPROX=ROM_hip_stance_R_PROX(sort_idx);
ROM_hip_swing_RsortedPROX=ROM_hip_swing_R_PROX(sort_idx);
ROM_hip_stance_LsortedPROX=ROM_hip_stance_L_PROX(sort_idx);
ROM_hip_swing_LsortedPROX=ROM_hip_swing_L_PROX(sort_idx);


ROM_hip_stance_RsortedDIST=ROM_hip_stance_R_DIST(sort_idx);
ROM_hip_swing_RsortedDIST=ROM_hip_swing_R_DIST(sort_idx);
ROM_hip_stance_LsortedDIST=ROM_hip_stance_L_DIST(sort_idx);
ROM_hip_swing_LsortedDIST=ROM_hip_swing_L_DIST(sort_idx);


figure;
x_positions = 1:length(sorted_frequencies);
bar_width = 0.3;
bar(x_positions - bar_width/2, ROM_knee_L_sortedPROX, bar_width, 'g');
hold on;
x_positions = 1:length(sorted_frequencies);
bar(x_positions + bar_width/2, ROM_knee_R_sortedPROX, bar_width, 'y');
xlabel('Frequencies');
ylabel('ROM_knee');
xticks(x_positions);
xticklabels(cellfun(@(x) [num2str(x) ' Hz'], num2cell(sorted_frequencies), 'UniformOutput', false));
title('R.O.M. KNEE PROX');
legend('Left knee', 'Right knee', 'Location', 'northwest');

figure;
x_positions = 1:length(sorted_frequencies);
bar_width = 0.3;
bar(x_positions - bar_width/2, ROM_knee_L_sortedDIST, bar_width, 'g');
hold on;
x_positions = 1:length(sorted_frequencies);
bar(x_positions + bar_width/2, ROM_knee_R_sortedDIST, bar_width, 'y');
xlabel('Frequencies');
ylabel('ROM_knee');
xticks(x_positions);
xticklabels(cellfun(@(x) [num2str(x) ' Hz'], num2cell(sorted_frequencies), 'UniformOutput', false));
title('R.O.M. KNEE DIST');
legend('Left knee', 'Right knee', 'Location', 'northwest');




figure;
x_positions = 1:length(sorted_frequencies);
bar_width = 0.3;
bar(x_positions - bar_width/2, ROM_hip_L_sortedPROX, bar_width, 'g');
hold on;
x_positions = 1:length(sorted_frequencies);
bar(x_positions + bar_width/2, ROM_hip_R_sortedPROX, bar_width, 'y');
xlabel('Frequencies');
ylabel('ROM_hip');
xticks(x_positions);
xticklabels(cellfun(@(x) [num2str(x) ' Hz'], num2cell(sorted_frequencies), 'UniformOutput', false));
title('R.O.M. HIP  PROX');
legend('Left hip', 'Right hip', 'Location', 'northwest');


figure;
x_positions = 1:length(sorted_frequencies);
bar_width = 0.3;
bar(x_positions - bar_width/2, ROM_hip_L_sortedDIST, bar_width, 'g');
hold on;
x_positions = 1:length(sorted_frequencies);
bar(x_positions + bar_width/2, ROM_hip_R_sortedDIST, bar_width, 'y');
xlabel('Frequencies');
ylabel('ROM_hip');
xticks(x_positions);
xticklabels(cellfun(@(x) [num2str(x) ' Hz'], num2cell(sorted_frequencies), 'UniformOutput', false));
title('R.O.M. HIP  DIST');
legend('Left hip', 'Right hip', 'Location', 'northwest');



%  %%  STATISTICAL ANALYSIS STANCE
% StatisticalAnalysis_stance(stance_duration_R);
% StatisticalAnalysis_stance(stance_duration_L);
% StatisticalAnalysis_stance(swing_duration_R);
% StatisticalAnalysis_stance(swing_duration_L);
% 
% StatisticalAnalysis_stance(gait_duration_R);
% StatisticalAnalysis_stance(gait_duration_L);
% 
% StatisticalAnalysis_stance(percentage_matrix_swingR);
% StatisticalAnalysis_stance(percentage_matrix_stanceR);
% %% STATISTICAL ANALYSIS ROM
% 
% StatisticalAnalysis_ROM(ROM_knee_all_L);
% StatisticalAnalysis_ROM(ROM_knee_all_R);
% StatisticalAnalysis_ROM(ROM_hip_all_L);
% StatisticalAnalysis_ROM(ROM_hip_all_R);
% 
% %% FULLGAIT FREQUENCY PLOTS
% % Nomi dei muscoli
% muscles_L = {'TA_L', 'VM_L', 'VL_L','RF_L', 'GM_L', 'ST_L', 'GASTRO_L', 'SO_L'};
% muscles_R = {'TA_R', 'VM_R', 'VL_R','RF_R', 'GM_R', 'ST_R', 'GASTRO_R', 'SO_R'};
% 
% % Frequenze corrispondenti (escludendo l'undicesima colonna)
% frequencies = [100, 10000, 15, 1000, 30, 50, 5, 5000, 75, 0];
% 
% % Rimuovere l'undicesima colonna dai dati
% max_values_L = max_values_L(:, 1:10);
% max_values_R = max_values_R(:, 1:10);
% 
% % Rimuovere l'undicesima colonna dalle frequenze
% frequencies = frequencies(1:10);
% 
% % Ordinare le frequenze e i corrispondenti valori di picco
% [sorted_frequencies, sort_idx] = sort(frequencies);
% 
% max_values_L_sorted = max_values_L(:, sort_idx);
% max_values_R_sorted = max_values_R(:, sort_idx);
% 
% % Parametri per il filtro di Savitzky-Golay
% sgolay_order = 1; % Ordine del polinomio
% sgolay_frame = 3; % Dimensione della finestra (più grande per un smoothing maggiore)
% 
% % Creazione delle figure per il plotting
% figure;
% 
% % Plot per i muscoli della gamba sinistra
% for i = 1:8
%     subplot(2, 4, i); % Creazione di subplot 2x4
%     
%     % Applica il filtro di Savitzky-Golay
%     smoothed_values_L = sgolayfilt(max_values_L_sorted(i, :), sgolay_order, sgolay_frame);
%     
%     semilogx(sorted_frequencies,smoothed_values_L, '-o', 'LineWidth', 2); % Plot logaritmico
%     hold on;
%     xline(100, '--r', 'LineWidth', 1.5); % Linea verticale tratteggiata rossa a 100 Hz
%     title(muscles_L{i});
%     xlabel('Frequency (Hz)');
%     ylabel('Peak Value');
%     grid on;
%     hold off;
% end
% 
% figure;
% 
% % Plot per i muscoli della gamba destra
% for i = 1:8
%     subplot(2, 4, i); % Creazione di subplot 2x4
%     
%     % Applica il filtro di Savitzky-Golay
%     smoothed_values_R = sgolayfilt(max_values_R_sorted(i, :), sgolay_order, sgolay_frame);
%     
%     semilogx(sorted_frequencies, smoothed_values_R, '-o', 'LineWidth', 2); % Plot logaritmico
%     hold on;
%     xline(100, '--r', 'LineWidth', 1.5); % Linea verticale tratteggiata rossa a 100 Hz
%     title(muscles_R{i});
%     xlabel('Frequency (Hz)');
%     ylabel('Peak Value');
%     grid on;
%     hold off;
% end
% 
% 
%%  STANCE & SWING


%PLOT SIGNALS OF RIGHT LEG
[max_values_Rstance]= Stance_Cal_subplot_R(emg_files, fs_IMU, fs_EMG);

% PLOT SIGNALS OF THE LEFT LEG
[max_values_Rswing]= Swing_Cal_subplot_R(emg_files, fs_IMU, fs_EMG);

%PLOT SIGNALS OF RIGHT LEG
[max_values_Lstance]= Stance_Cal_subplot_L(emg_files, fs_IMU, fs_EMG);

% PLOT SIGNALS OF THE LEFT LEG
[max_values_Lswing]= Swing_Cal_subplot_L(emg_files, fs_IMU, fs_EMG);






%% STANCE/SWING RIGHT FREQUENCY PLOTS

 muscles_R = {'GM_R', 'GASTRO_R','TA_R', 'ST_R', 'VM_R', 'VL_R', 'RF_R'};
 % Frequencies of the corrisponding value
 frequencies = [100,15,30,50,5,65,75,0];
 
% Loop su ogni riga della matrice
for i = 1:size(max_values_Rstance, 1)
    % Controlla se almeno un elemento della riga è superiore a 100
    if any(max_values_Rstance(i, :) > 100)
        % Trova il massimo valore nella riga
        max_value = max(max_values_Rstance(i, :));
        
        % Normalizza la riga dividendo ogni elemento per il massimo
        max_values_Rstance(i, :) = (max_values_Rstance(i, :) / max_value)*100;
    end
end


% Loop su ogni riga della matrice
for i = 1:size(max_values_Rswing, 1)
    % Controlla se almeno un elemento della riga è superiore a 100
    if any(max_values_Rswing(i, :) > 100)
        % Trova il massimo valore nella riga
        max_value = max(max_values_Rswing(i, :));
        
        % Normalizza la riga dividendo ogni elemento per il massimo
        max_values_Rswing(i, :) = (max_values_Rswing(i, :) / max_value)*100;
    end
end


% Delete the last column
 max_values_RstanceDIST = max_values_Rstance(:,[1:7,15]);
 max_values_RswingDIST = max_values_Rswing(:,[1:7,15]);

 max_values_RstancePROX = max_values_Rstance(:,8:15);
 max_values_RswingPROX = max_values_Rswing(:,8:15);
 
% Sort the frequencies and relative peak values
 [sorted_frequencies, sort_idx] = sort(frequencies);
% 
 max_values_Rstance_sortedPROX = max_values_RstancePROX(:, sort_idx);
 max_values_Rswing_sortedPROX = max_values_RswingPROX(:, sort_idx);

 max_values_Rstance_sortedDIST = max_values_RstanceDIST(:, sort_idx);
 max_values_Rswing_sortedDIST = max_values_RswingDIST(:, sort_idx);
 
%  Savitzky-Golay filter parameters
sgolay_order = 1; % Polynomial order
sgolay_frame = 3; % Window lenght (bigger, bigger smooth action)


figure;
for i = 1:7
    subplot(2, 4, i); 
    
    % Apply Savitzky-Golay filter
    smoothed_values_stance = sgolayfilt(max_values_Rstance_sortedPROX(i, :), sgolay_order, sgolay_frame);
    %semilogx(sorted_frequencies,smoothed_values_stance, '-o', 'LineWidth', 2); % Plot logaritmico
    plot(sorted_frequencies,smoothed_values_stance, '-o', 'LineWidth', 2); % Plot logaritmico
    hold on;
    title('Stance phase PROX', muscles_R{i});
    xlabel('Frequency (Hz)');
    ylabel('Peak Value');
    grid on;
    hold off;
end

figure;
for i = 1:7
    subplot(2, 4, i); 
    
    %Apply Savitzky-Golay filter
    smoothed_values_swing = sgolayfilt(max_values_Rswing_sortedPROX(i, :), sgolay_order, sgolay_frame);
    %semilogx(sorted_frequencies, smoothed_values_swing, '-o', 'LineWidth', 2); % Plot logaritmico
    plot(sorted_frequencies, smoothed_values_swing, '-o', 'LineWidth', 2); % Plot logaritmico
    hold on;
%     xline(100, '--r', 'LineWidth', 1.5); 
    title('Swing phase PROX',muscles_R{i});
    xlabel('Frequency (Hz)');
    ylabel('Peak Value');
    grid on;
    hold off;
end

figure;
for i = 1:7
    subplot(2, 4, i); 
    
    % Apply Savitzky-Golay filter
    smoothed_values_stance = sgolayfilt(max_values_Rstance_sortedDIST(i, :), sgolay_order, sgolay_frame);
    %semilogx(sorted_frequencies,smoothed_values_stance, '-o', 'LineWidth', 2); % Plot logaritmico
    plot(sorted_frequencies,smoothed_values_stance, '-o', 'LineWidth', 2); % Plot logaritmic
    hold on;
%     xline(100, '--r', 'LineWidth', 1.5); % Linea verticale tratteggiata rossa a 100 Hz
    title('Stance phase DIST', muscles_R{i});
    xlabel('Frequency (Hz)');
    ylabel('Peak Value');
    grid on;
    hold off;
end

figure;
for i = 1:7
    subplot(2, 4, i); 
    
    %Apply Savitzky-Golay filter
    smoothed_values_swing = sgolayfilt(max_values_Rswing_sortedDIST(i, :), sgolay_order, sgolay_frame);
    %semilogx(sorted_frequencies, smoothed_values_swing, '-o', 'LineWidth', 2); % Plot logaritmico
    plot(sorted_frequencies, smoothed_values_swing, '-o', 'LineWidth', 2); % Plot logaritmico
    hold on;
    %xline(100, '--r', 'LineWidth', 1.5); 
    title('Swing phase DIST',muscles_R{i});
    xlabel('Frequency (Hz)');
    ylabel('Peak Value');
    grid on;
    hold off;
end

%% STANCE/SWING LEFT FREQUENCY PLOTS

muscles_L = {'GASTRO_L', 'TA_L','ST_L' 'VM_L', 'VL_L', 'RF_L', 'GM_L'};

% Frequencies of the corrisponding value
frequencies = [100,15,30,50,5,65,75,0];


% Loop su ogni riga della matrice
for i = 1:size(max_values_Lstance, 1)
    % Controlla se almeno un elemento della riga è superiore a 100
    if any(max_values_Lstance(i, :) > 100)
        % Trova il massimo valore nella riga
        max_value = max(max_values_Lstance(i, :));
        
        % Normalizza la riga dividendo ogni elemento per il massimo
        max_values_Lstance(i, :) = (max_values_Lstance(i, :) / max_value)*100;
    end
end


% Loop su ogni riga della matrice
for i = 1:size(max_values_Lswing, 1)
    % Controlla se almeno un elemento della riga è superiore a 100
    if any(max_values_Lswing(i, :) > 100)
        % Trova il massimo valore nella riga
        max_value = max(max_values_Lswing(i, :));
        
        % Normalizza la riga dividendo ogni elemento per il massimo
        max_values_Lswing(i, :) = (max_values_Lswing(i, :) / max_value)*100;
    end
end


% Delete the last column
 max_values_LstancePROX = max_values_Lstance(:, [1:7,15]);
 max_values_LswingPROX = max_values_Lswing(:, [1:7,15]);

 max_values_LstanceDIST = max_values_Lstance(:, 8:15);
 max_values_LswingDIST = max_values_Lswing(:, 8:15);

 % Sort the frequencies and relative peak values
 [sorted_frequencies, sort_idx] = sort(frequencies);
 
 max_values_Lstance_sortedPROX = max_values_LstancePROX(:, sort_idx);
 max_values_Lswing_sortedPROX = max_values_LswingPROX(:, sort_idx);
 max_values_Lstance_sortedDIST = max_values_LstanceDIST(:, sort_idx);
 max_values_Lswing_sortedDIST = max_values_LswingDIST(:, sort_idx);

 
% Savitzky-Golay Filter parameters
sgolay_order = 1; % Polynomial order
sgolay_frame = 3; % Window length (bigger, bigger smooth action)

figure;
for i = 1:7
    subplot(2, 4, i); 
    
    % Apply Savitzky-Golay filter
    smoothed_values_stance = sgolayfilt(max_values_Lstance_sortedPROX(i, :), sgolay_order, sgolay_frame);
    plot(sorted_frequencies,smoothed_values_stance, '-o', 'LineWidth', 2); % Plot logaritmico
    hold on;
    %xline(100, '--r', 'LineWidth', 1.5); % Linea verticale tratteggiata rossa a 100 Hz
    title('Stance phase PROX', muscles_L{i});
    xlabel('Frequency (Hz)');
    ylabel('Peak Value');
    grid on;
    hold off;
end

figure;
for i = 1:7
    subplot(2, 4, i);
    
    % Apply Savitzky-Golay filter
    smoothed_values_swing = sgolayfilt(max_values_Lswing_sortedPROX(i, :), sgolay_order, sgolay_frame);
    plot(sorted_frequencies, smoothed_values_swing, '-o', 'LineWidth', 2); % Plot logaritmico
    hold on;
    %xline(100, '--r', 'LineWidth', 1.5); 
    title('Swing phase PROX',muscles_L{i});
    xlabel('Frequency (Hz)');
    ylabel('Peak Value');
    grid on;
    hold off;
end


figure;
for i = 1:7
    subplot(2, 4, i); 
    
    % Apply Savitzky-Golay filter
    smoothed_values_stance = sgolayfilt(max_values_Lstance_sortedDIST(i, :), sgolay_order, sgolay_frame);
    plot(sorted_frequencies,smoothed_values_stance, '-o', 'LineWidth', 2); % Plot logaritmico
    hold on;
    %xline(100, '--r', 'LineWidth', 1.5); % Linea verticale tratteggiata rossa a 100 Hz
    title('Stance phase DIST', muscles_L{i});
    xlabel('Frequency (Hz)');
    ylabel('Peak Value');
    grid on;
    hold off;
end

figure;
for i = 1:7
    subplot(2, 4, i);
    
    % Apply Savitzky-Golay filter
    smoothed_values_swing = sgolayfilt(max_values_Lswing_sortedDIST(i, :), sgolay_order, sgolay_frame);
    plot(sorted_frequencies, smoothed_values_swing, '-o', 'LineWidth', 2); % Plot logaritmico
    hold on;
    %xline(100, '--r', 'LineWidth', 1.5); 
    title('Swing phase DIST',muscles_L{i});
    xlabel('Frequency (Hz)');
    ylabel('Peak Value');
    grid on;
    hold off;
end





%% LEFT + RIGHT VALUES TOGETHER
 
 gait_duration_TOT=[gait_duration_R,gait_duration_L]; 

% Specifica il nome del file con cui salvare la matrice
% fileName_gait =  'GAIT_DURATION14';  % Sostituisci 'nome_matrice' con il nome desiderato
% fileName_stance= 'STANCE_DURATION14';
% fileName_swing = 'SWING_DURATION14';
% FileName_emgRstance='EMG_RSTANCE14';
% FileName_emgLstance='EMG_LSTANCE14';
% FileName_emgRswing='EMG_RSWING14';
% FileName_emgLswing='EMG_LSWING14';


% Combina il base path con il nome del file
% fullPath_gait = fullfile('C:\Users\39320\Desktop\EPFL\Thesis code\tSCS_Continuos\NEW_PROTOCOL\MATRIX\GAIT_DURATION', [fileName_gait, '.mat']);
% fullPath_stance = fullfile('C:\Users\39320\Desktop\EPFL\Thesis code\tSCS_Continuos\NEW_PROTOCOL\MATRIX\STANCE_DURATION', [fileName_stance, '.mat']);
% fullPath_swing = fullfile('C:\Users\39320\Desktop\EPFL\Thesis code\tSCS_Continuos\NEW_PROTOCOL\MATRIX\SWING_DURATION', [fileName_swing, '.mat']);
% fullPath_emgRstance = fullfile('C:\Users\39320\Desktop\EPFL\Thesis code\tSCS_Continuos\NEW_PROTOCOL\MATRIX\EMG', [FileName_emgRstance, '.mat']);
% fullPath_emgLstance = fullfile('C:\Users\39320\Desktop\EPFL\Thesis code\tSCS_Continuos\NEW_PROTOCOL\MATRIX\EMG', [FileName_emgLstance, '.mat']);
% fullPath_emgRswing = fullfile('C:\Users\39320\Desktop\EPFL\Thesis code\tSCS_Continuos\NEW_PROTOCOL\MATRIX\EMG', [FileName_emgRswing, '.mat']);
% fullPath_emgLswing = fullfile('C:\Users\39320\Desktop\EPFL\Thesis code\tSCS_Continuos\NEW_PROTOCOL\MATRIX\EMG', [FileName_emgLswing, '.mat']);

% Salva la matrice usando il percorso completo
% save(fullPath_gait, 'gait_duration_TOT');
% save(fullPath_emgRstance, 'max_values_Rstance');
% save(fullPath_emgLstance, 'max_values_Lstance');
% save(fullPath_emgRswing, 'max_values_Rswing');
% save(fullPath_emgLswing, 'max_values_Lswing');

%Preallocation
numRows = size(gait_duration_TOT, 1);
medie = zeros(numRows, 1);

for i = 1:numRows
    % Select only the elements =! 0
    elementiValidi = gait_duration_TOT(i, gait_duration_TOT(i, :) ~= 0);
    
    % Mean calculation
    if ~isempty(elementiValidi)
        medie(i) = mean(elementiValidi);
    else
        medie(i) = NaN; 
    end
end

% Visualize mean
disp(medie);
gait_meanPROX=medie(8:15);
gait_meanDIST=medie([1:7,15]);


gait_meanPROXsorted=gait_meanPROX(sort_idx)';
gait_meanDISTsorted=gait_meanDIST(sort_idx)';


stance_duration_TOT=[stance_duration_R,stance_duration_L]; 
% Salva la matrice usando il percorso completo
% save(fullPath_stance, 'stance_duration_TOT');



%Preallocation
numRows = size(stance_duration_TOT, 1);
medie_stance = zeros(numRows, 1);

for i = 1:numRows
    % Select only the elements =! 0
    elementiValidi_stance = stance_duration_TOT(i, stance_duration_TOT(i, :) ~= 0);
    
    % Mean calculation
    if ~isempty(elementiValidi_stance)
        medie_stance(i) = mean(elementiValidi_stance);
    else
        medie_stance(i) = NaN; 
    end
end

% Visualize mean
disp(medie_stance);
stance_meanPROX=medie_stance(8:15);
stance_meanDIST=medie_stance([1:7,15]);


stance_meanPROXsorted=stance_meanPROX(sort_idx)';
stance_meanDISTsorted=stance_meanDIST(sort_idx)';


swing_duration_TOT=[swing_duration_R,swing_duration_L]; 
% save(fullPath_swing, 'swing_duration_TOT');


%Preallocation
numRows = size(swing_duration_TOT, 1);
medie_swing = zeros(numRows, 1);

for i = 1:numRows
    % Select only the elements =! 0
    elementiValidi_swing = swing_duration_TOT(i, swing_duration_TOT(i, :) ~= 0);
    
    % Mean calculation
    if ~isempty(elementiValidi_swing)
        medie_swing(i) = mean(elementiValidi_swing);
    else
        medie_swing(i) = NaN; 
    end
end

% Visualize mean
disp(medie_swing);
swing_meanPROX=medie_swing(8:15);
swing_meanDIST=medie_swing([1:7,15]);


swing_meanPROXsorted=swing_meanPROX(sort_idx)';
swing_meanDISTsorted=swing_meanDIST(sort_idx)';




%  for i = 1:11
%     ROM_hip_all_TOT{i} = [ROM_hip_all_L{i}; ROM_hip_all_R{i}]; 
%  end
% 
%   for i = 1:11
%     ROM_knee_all_TOT{i} = [ROM_knee_all_L{i}; ROM_knee_all_R{i}]; 
%   end
% 
% for i = 1:11
%     HIPmean_TOT(i) = mean(ROM_hip_all_TOT{i});
% end
% disp(HIPmean_TOT);
% 
% for i = 1:11
%     KNEEmean_TOT(i) = mean(ROM_knee_all_TOT{i});
% end
% disp(KNEEmean_TOT);

%% HIP SWING

 FileName_HIPSWING='HIPSWING3';
 FileName_HIPSTANCE='HIPSTANCE3';
 FileName_KNEESWING='KNEESWING3';
 FileName_KNEESTANCE='KNEESTANCE3';


 fullPath_HIPSWING = fullfile('C:\Users\39320\Desktop\EPFL\Thesis code\tSCS_Continuos\NEW_PROTOCOL\MATRIX\HIP_SWING', [FileName_HIPSWING, '.mat']);
 fullPath_HIPSTANCE = fullfile('C:\Users\39320\Desktop\EPFL\Thesis code\tSCS_Continuos\NEW_PROTOCOL\MATRIX\HIP_STANCE', [FileName_HIPSTANCE, '.mat']);
 fullPath_KNEESWING = fullfile('C:\Users\39320\Desktop\EPFL\Thesis code\tSCS_Continuos\NEW_PROTOCOL\MATRIX\KNEE_SWING', [FileName_KNEESWING, '.mat']);
 fullPath_KNEESTANCE = fullfile('C:\Users\39320\Desktop\EPFL\Thesis code\tSCS_Continuos\NEW_PROTOCOL\MATRIX\KNEE_STANCE', [FileName_KNEESTANCE, '.mat']);




% Trova la lunghezza massima dei vettori nelle celle
maxLengthHIPSW_L = max(cellfun(@length, ROM_hip_swing_all_L));
maxLengthHIPSW_R = max(cellfun(@length, ROM_hip_swing_all_R));

% Inizializza la matrice con NaN o zeri (a seconda delle tue esigenze)
matrixHIPSW_L = zeros(15, maxLengthHIPSW_L); % Puoi usare 0 se preferisci

% Riempie la matrice
for i = 1:15
    vector = ROM_hip_swing_all_L{i};
    matrixHIPSW_L(i, 1:length(vector)) = vector;
end


% Inizializza la matrice con NaN o zeri (a seconda delle tue esigenze)
matrixHIPSW_R = zeros(15, maxLengthHIPSW_R); % Puoi usare 0 se preferisci

% Riempie la matrice
for i = 1:15
    vector = ROM_hip_swing_all_R{i};
    matrixHIPSW_R(i, 1:length(vector)) = vector;
end

matrixHIPSW_TOT=[matrixHIPSW_R,matrixHIPSW_L];
 save(fullPath_HIPSWING, 'matrixHIPSW_TOT');


% HIP STANCE  
% Trova la lunghezza massima dei vettori nelle celle
maxLengthHIPST_L = max(cellfun(@length, ROM_hip_stance_all_L));
maxLengthHIPST_R = max(cellfun(@length, ROM_hip_stance_all_R));

% Inizializza la matrice con NaN o zeri (a seconda delle tue esigenze)
matrixHIPST_L = zeros(15, maxLengthHIPST_L); % Puoi usare 0 se preferisci

% Riempie la matrice
for i = 1:15
    vector = ROM_hip_stance_all_L{i};
    matrixHIPST_L(i, 1:length(vector)) = vector;
end


% Inizializza la matrice con NaN o zeri (a seconda delle tue esigenze)
matrixHIPST_R = zeros(15, maxLengthHIPST_R); % Puoi usare 0 se preferisci

% Riempie la matrice
for i = 1:15
    vector = ROM_hip_stance_all_R{i};
    matrixHIPST_R(i, 1:length(vector)) = vector;
end

matrixHIPST_TOT=[matrixHIPST_R,matrixHIPST_L];
save(fullPath_HIPSTANCE, 'matrixHIPST_TOT');

% KNEE SWING


% Trova la lunghezza massima dei vettori nelle celle
maxLengthKNEESW_L = max(cellfun(@length, ROM_knee_swing_all_L));
maxLengthKNEESW_R = max(cellfun(@length, ROM_knee_swing_all_R));

% Inizializza la matrice con NaN o zeri (a seconda delle tue esigenze)
matrixKNEESW_L = zeros(15, maxLengthKNEESW_L); % Puoi usare 0 se preferisci

% Riempie la matrice
for i = 1:15
    vector = ROM_knee_swing_all_L{i};
    matrixKNEESW_L(i, 1:length(vector)) = vector;
end


% Inizializza la matrice con NaN o zeri (a seconda delle tue esigenze)
matrixKNEESW_R = zeros(15, maxLengthKNEESW_R); % Puoi usare 0 se preferisci

% Riempie la matrice
for i = 1:15
    vector = ROM_knee_swing_all_R{i};
    matrixKNEESW_R(i, 1:length(vector)) = vector;
end

matrixKNEESW_TOT=[matrixKNEESW_R,matrixKNEESW_L];
save(fullPath_KNEESWING, 'matrixKNEESW_TOT');


% KNEE STANCE  
% Trova la lunghezza massima dei vettori nelle celle
maxLengthKNEEST_L = max(cellfun(@length, ROM_knee_stance_all_L));
maxLengthKNEEST_R = max(cellfun(@length, ROM_knee_stance_all_R));

% Inizializza la matrice con NaN o zeri (a seconda delle tue esigenze)
matrixKNEEST_L = zeros(15, maxLengthKNEEST_L); % Puoi usare 0 se preferisci

% Riempie la matrice
for i = 1:15
    vector = ROM_knee_stance_all_L{i};
    matrixKNEEST_L(i, 1:length(vector)) = vector;
end


% Inizializza la matrice con NaN o zeri (a seconda delle tue esigenze)
matrixKNEEST_R = zeros(15, maxLengthKNEEST_R); % Puoi usare 0 se preferisci

% Riempie la matrice
for i = 1:15
    vector = ROM_knee_stance_all_R{i};
    matrixKNEEST_R(i, 1:length(vector)) = vector;
end

matrixKNEEST_TOT=[matrixKNEEST_R,matrixKNEEST_L];
save(fullPath_KNEESTANCE, 'matrixKNEEST_TOT');