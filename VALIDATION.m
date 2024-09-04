clear all
close all
clc


load  'nello xsens validation'\XSENSvalidation\nello_1.mat

%% RIGHT

% Frequenza di campionamento
fs_IMU = 148.1481475830078; % Frequenza di campionamento in Hz

% Calcola i tempi in secondi corrispondenti ai campioni di dati
tempi_secondi = (0:length(Data(2,:))-1) / fs_IMU;

% Definisci la lunghezza della finestra per la media mobile
    window_size = 12;  % Numero di campioni

    % Applica il filtro a media mobile
    data_filtered = movmean(Data(2,:), window_size);
    

% Trova i picchi di accelerazione verticale (indicativi del contatto del tallone)
[peaks, peak_locs] = findpeaks(Data(2,:), 'MinPeakHeight', 0.25, 'MinPeakDistance', 0.8 * 148.15); % Imposta una distanza minima di 1.5 secondi

% Trova i minimi locali di accelerazione verticale (indicativi del sollevamento del tallone)
[troughs, trough_locs] = findpeaks(-Data(2,:), 'MinPeakHeight', 0.25, 'MinPeakDistance', 0.8 * 148.15); % Ricerca i minimi locali invertendo il segnale e trovando i picchi

num_contatti_tallone = length(peaks)
num_sollevamenti_tallone = length(troughs)

% Plot dei dati e dei punti di riferimento identificati
figure;
plot(tempi_secondi, Data(2,:)); hold on;
plot(tempi_secondi(peak_locs), peaks, 'ro', 'MarkerSize', 10); % Picchi (contatto del tallone)
plot(tempi_secondi(trough_locs), -troughs, 'go', 'MarkerSize', 10); % Minimi locali (sollevamento del tallone)
xlabel('Tempo (s)');
ylabel('Accelerazione verticale');
legend('Accelerazione verticale', 'Contatto del tallone', 'Sollevamento del tallone');
title('Identificazione dei punti di riferimento del ciclo di cammino');

% Limita l'asse x fino a 350 secondi
xlim([0, 350]);

% Imposta i limiti dell'asse y in base ai valori massimi e minimi dei dati
% Calcola il valore massimo e minimo dei dati di accelerazione verticale
max_accelerazione = max(abs(Data(2,:)));
min_accelerazione = -max_accelerazione;

% Calcola i limiti dell'asse y aggiungendo un margine percentuale
margine_percentuale = 0.1; % 10% di margine
margine_superiore = max_accelerazione * (1 + margine_percentuale);
margine_inferiore = min_accelerazione * (1 - margine_percentuale);
ylim([margine_inferiore, margine_superiore]);


% % Individua gli istanti di tempo corrispondenti ai contatti del tallone e ai sollevamenti
tempi_contatti_tallone_R = tempi_secondi(peak_locs);
tempi_sollevamenti_tallone_R = tempi_secondi(trough_locs);


    figure;
    ROM_knee_R=[];
    ROM_hip_R=[];
    ROM_hip_all_R =[];
    ROM_knee_all_R = []; 

        
       f=length(Data(5,:));
       gyro_tibia_signal = Data(5,1:(f/8.5));

       h=length(Data(12,:));
       gyro_cosc_signal = Data(12,1:(h/8.5));

       l=length(Data(19,:));
       gyro_bacino_signal= Data(19,1:(l/8.5));

       num_cycles = length(tempi_contatti_tallone_R) - 1;

       max_cycle_length = 0;
for j = 1:num_cycles
    start_time = tempi_contatti_tallone_R(j);
    end_time = tempi_contatti_tallone_R(j+1);
    cycle_length = end_time - start_time;
    if cycle_length > max_cycle_length
        max_cycle_length = cycle_length;
    end
end
theta_anca_all_cycles = [];
theta_ginocchio_all_cycles = [];

% Loop su ogni ciclo di camminata
for h = 1:num_cycles
    start_time = tempi_contatti_tallone_R(h);
    end_time = tempi_contatti_tallone_R(h+1);

    % Conversione dell'intervallo di tempo in indici di campionamento
    start_index = round(start_time * fs_IMU);
    end_index = round(end_time * fs_IMU);

    % Estrazione del segnale nell'intervallo specificato
    gyro_tibia_segment = gyro_tibia_signal(start_index:end_index);
    gyro_cosc_segment = gyro_cosc_signal(start_index:end_index);
    gyro_bacino_segment= gyro_bacino_signal(start_index:end_index);

    % Calcolo di dt
    dt = 1 / fs_IMU;

    % Integrazione dei dati del giroscopio filtrati per ottenere gli angoli
    theta_tibia = cumtrapz((0:length(gyro_tibia_segment)-1)*dt, gyro_tibia_segment); % Angolo della tibia
    theta_cosc = cumtrapz((0:length(gyro_cosc_segment)-1)*dt, gyro_cosc_segment); % Angolo della coscia
    theta_bacino = cumtrapz((0:length(gyro_bacino_segment)-1)*dt, gyro_bacino_segment); % Angolo della coscia

    % Calcolo dell'angolo del ginocchio come differenza tra angoli di tibia e coscia
    theta_ginocchio = theta_cosc - theta_tibia;
    theta_anca=  theta_bacino + theta_cosc;

    % Normalizzazione del tempo
    normalized_time_knee = linspace(0, 1, length(gyro_tibia_segment));

    % Interpolazione dei segnali per avere la stessa lunghezza temporale
    normalized_theta_ginocchio = interp1(normalized_time_knee, theta_ginocchio, linspace(0, 1, round(max_cycle_length * fs_IMU)));

    % Memorizzazione degli angoli normalizzati del ginocchio
    theta_ginocchio_all_cycles = [theta_ginocchio_all_cycles; normalized_theta_ginocchio];

subplot(2,1,1);
plot(linspace(0, 100, size(theta_ginocchio_all_cycles, 2)), theta_ginocchio_all_cycles, 'LineWidth', 2);
title('knee angle mean Right');


     % Normalizzazione del tempo
    normalized_time_hip = linspace(0, 1, length(gyro_bacino_segment));

    % Interpolazione dei segnali per avere la stessa lunghezza temporale
    normalized_theta_anca = interp1(normalized_time_hip, theta_anca, linspace(0, 1, round(max_cycle_length * fs_IMU)));

    % Memorizzazione degli angoli normalizzati dell'anca
    theta_anca_all_cycles = [theta_anca_all_cycles; normalized_theta_anca];
  subplot(2,1,2);
plot(linspace(0, 100, size(theta_anca_all_cycles, 2)), theta_anca_all_cycles, 'LineWidth', 2);
title('Hip angle mean Right');
end

figure;
mean_theta_anca = mean(theta_anca_all_cycles, 1);
mean_theta_ginocchio = mean(theta_ginocchio_all_cycles, 1);

ROM_MEAN_RIGHT_HIP=(max(mean_theta_anca)-min(mean_theta_anca));
ROM_MEAN_RIGHT_KNEE=(max(mean_theta_ginocchio)-min(mean_theta_ginocchio));

hold on;
% Plot della media degli angoli dell'anca
subplot(2,1,1);
plot(linspace(0, 100, size(mean_theta_anca, 2)), mean_theta_anca, 'LineWidth', 2);
title('Hip angle mean Right');
xlabel('Gait cycle % ');

hold on;
subplot(2,1,2);
plot(linspace(0, 100, size(mean_theta_ginocchio, 2)), mean_theta_ginocchio, 'LineWidth', 2);
title('Knee angle mean Right');
xlabel('Gait cycle %');


% Calcolo del Range of Motion (ROM) per l'anca per ogni ciclo
ROM_hip_all_R_i = max(theta_anca_all_cycles, [], 2) - min(theta_anca_all_cycles, [], 2);
ROM_hip_R = mean(ROM_hip_all_R_i);
ROM_hip_all_R = ROM_hip_all_R_i; % Store ROM_hip_all_L data in cell


% Calcolo del Range of Motion (ROM) per il ginocchio per ogni ciclo
ROM_knee_all_R_i = max(theta_ginocchio_all_cycles, [], 2) - min(theta_ginocchio_all_cycles, [], 2);
ROM_knee_R = mean(ROM_knee_all_R_i);
ROM_knee_all_R = ROM_knee_all_R_i;


%% LEFT
% Frequenza di campionamento
fs_IMU = 148.1481475830078; % Frequenza di campionamento in Hz

% Calcola i tempi in secondi corrispondenti ai campioni di dati
tempi_secondi = (0:length(Data(23,:))-1) / fs_IMU;

% Definisci la lunghezza della finestra per la media mobile
    window_size = 12;  % Numero di campioni

    % Applica il filtro a media mobile
    data_filtered = movmean(Data(23,:), window_size);
    

% Trova i picchi di accelerazione verticale (indicativi del contatto del tallone)
[peaks, peak_locs] = findpeaks(Data(23,:), 'MinPeakHeight', 0.25, 'MinPeakDistance', 0.8 * 148.15); % Imposta una distanza minima di 1.5 secondi

% Trova i minimi locali di accelerazione verticale (indicativi del sollevamento del tallone)
[troughs, trough_locs] = findpeaks(-Data(23,:), 'MinPeakHeight', 0.25, 'MinPeakDistance', 0.8 * 148.15); % Ricerca i minimi locali invertendo il segnale e trovando i picchi

num_contatti_tallone = length(peaks)
num_sollevamenti_tallone = length(troughs)

% Plot dei dati e dei punti di riferimento identificati
figure;
plot(tempi_secondi, Data(23,:)); hold on;
plot(tempi_secondi(peak_locs), peaks, 'ro', 'MarkerSize', 10); % Picchi (contatto del tallone)
plot(tempi_secondi(trough_locs), -troughs, 'go', 'MarkerSize', 10); % Minimi locali (sollevamento del tallone)
xlabel('Tempo (s)');
ylabel('Accelerazione verticale');
legend('Accelerazione verticale', 'Contatto del tallone', 'Sollevamento del tallone');
title('Identificazione dei punti di riferimento del ciclo di cammino');

% Limita l'asse x fino a 350 secondi
xlim([0, 350]);

% Imposta i limiti dell'asse y in base ai valori massimi e minimi dei dati
% Calcola il valore massimo e minimo dei dati di accelerazione verticale
max_accelerazione = max(abs(Data(23,:)));
min_accelerazione = -max_accelerazione;

% Calcola i limiti dell'asse y aggiungendo un margine percentuale
margine_percentuale = 0.1; % 10% di margine
margine_superiore = max_accelerazione * (1 + margine_percentuale);
margine_inferiore = min_accelerazione * (1 - margine_percentuale);
ylim([margine_inferiore, margine_superiore]);


% % Individua gli istanti di tempo corrispondenti ai contatti del tallone e ai sollevamenti
tempi_contatti_tallone_R = tempi_secondi(peak_locs);
tempi_sollevamenti_tallone_R = tempi_secondi(trough_locs);


    figure;
    ROM_knee_R=[];
    ROM_hip_R=[];
    ROM_hip_all_R =[];
    ROM_knee_all_R = []; 

        
       f=length(Data(26,:));
       gyro_tibia_signal = Data(26,1:(f/8.5));

       h=length(Data(33,:));
       gyro_cosc_signal = Data(33,1:(h/8.5));

       l=length(Data(40,:));
       gyro_bacino_signal= Data(40,1:(l/8.5));

       num_cycles = length(tempi_contatti_tallone_R) - 1;

       max_cycle_length = 0;
for j = 1:num_cycles
    start_time = tempi_contatti_tallone_R(j);
    end_time = tempi_contatti_tallone_R(j+1);
    cycle_length = end_time - start_time;
    if cycle_length > max_cycle_length
        max_cycle_length = cycle_length;
    end
end
theta_anca_all_cycles = [];
theta_ginocchio_all_cycles = [];

% Loop su ogni ciclo di camminata
for h = 1:num_cycles
    start_time = tempi_contatti_tallone_R(h);
    end_time = tempi_contatti_tallone_R(h+1);

    % Conversione dell'intervallo di tempo in indici di campionamento
    start_index = round(start_time * fs_IMU);
    end_index = round(end_time * fs_IMU);

    % Estrazione del segnale nell'intervallo specificato
    gyro_tibia_segment = gyro_tibia_signal(start_index:end_index);
    gyro_cosc_segment = gyro_cosc_signal(start_index:end_index);
    gyro_bacino_segment= gyro_bacino_signal(start_index:end_index);

    % Calcolo di dt
    dt = 1 / fs_IMU;

    % Integrazione dei dati del giroscopio filtrati per ottenere gli angoli
    theta_tibia = cumtrapz((0:length(gyro_tibia_segment)-1)*dt, gyro_tibia_segment); % Angolo della tibia
    theta_cosc = cumtrapz((0:length(gyro_cosc_segment)-1)*dt, gyro_cosc_segment); % Angolo della coscia
    theta_bacino = cumtrapz((0:length(gyro_bacino_segment)-1)*dt, gyro_bacino_segment); % Angolo della coscia

    % Calcolo dell'angolo del ginocchio come differenza tra angoli di tibia e coscia
    theta_ginocchio = theta_cosc - theta_tibia;
    theta_anca=  theta_bacino + theta_cosc;

    % Normalizzazione del tempo
    normalized_time_knee = linspace(0, 1, length(gyro_tibia_segment));

    % Interpolazione dei segnali per avere la stessa lunghezza temporale
    normalized_theta_ginocchio = interp1(normalized_time_knee, theta_ginocchio, linspace(0, 1, round(max_cycle_length * fs_IMU)));

    % Memorizzazione degli angoli normalizzati del ginocchio
    theta_ginocchio_all_cycles = [theta_ginocchio_all_cycles; normalized_theta_ginocchio];

subplot(2,1,1);
plot(linspace(0, 100, size(theta_ginocchio_all_cycles, 2)), theta_ginocchio_all_cycles, 'LineWidth', 2);
title('knee angle mean Right');


     % Normalizzazione del tempo
    normalized_time_hip = linspace(0, 1, length(gyro_bacino_segment));

    % Interpolazione dei segnali per avere la stessa lunghezza temporale
    normalized_theta_anca = interp1(normalized_time_hip, theta_anca, linspace(0, 1, round(max_cycle_length * fs_IMU)));

    % Memorizzazione degli angoli normalizzati dell'anca
    theta_anca_all_cycles = [theta_anca_all_cycles; normalized_theta_anca];
  subplot(2,1,2);
plot(linspace(0, 100, size(theta_anca_all_cycles, 2)), theta_anca_all_cycles, 'LineWidth', 2);
title('Hip angle mean Right');
end

figure;
mean_theta_anca = mean(theta_anca_all_cycles, 1);
mean_theta_ginocchio = mean(theta_ginocchio_all_cycles, 1);


ROM_MEAN_LEFT_HIP=(max(mean_theta_anca)-min(mean_theta_anca));
ROM_MEAN_LEFT_KNEE=(max(mean_theta_ginocchio)-min(mean_theta_ginocchio));
hold on;
% Plot della media degli angoli dell'anca
subplot(2,1,1);
plot(linspace(0, 100, size(mean_theta_anca, 2)), mean_theta_anca, 'LineWidth', 2);
title('Hip angle mean Left');
xlabel('Gait cycle % ');

hold on;
subplot(2,1,2);
plot(linspace(0, 100, size(mean_theta_ginocchio, 2)), mean_theta_ginocchio, 'LineWidth', 2);
title('Knee angle mean Left');
xlabel('Gait cycle %');


% Calcolo del Range of Motion (ROM) per l'anca per ogni ciclo
ROM_hip_all_R_i = max(theta_anca_all_cycles, [], 2) - min(theta_anca_all_cycles, [], 2);
ROM_hip_R = mean(ROM_hip_all_R_i);
ROM_hip_all_R = ROM_hip_all_R_i; % Store ROM_hip_all_L data in cell


% Calcolo del Range of Motion (ROM) per il ginocchio per ogni ciclo
ROM_knee_all_R_i = max(theta_ginocchio_all_cycles, [], 2) - min(theta_ginocchio_all_cycles, [], 2);
ROM_knee_R = mean(ROM_knee_all_R_i);
ROM_knee_all_R = ROM_knee_all_R_i;




