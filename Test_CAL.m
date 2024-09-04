
%% Estrazione ciclo di camminata da dati IMU

clear all
close all
clc


load DATA\SUB1\SUB1_Withoutstim1.mat;

% Frequenza di campionamento
fs_IMU = 148.1481475830078; % Frequenza di campionamento in Hz

% Calcola i tempi in secondi corrispondenti ai campioni di dati
tempi_secondi = (0:length(Data(52,:))-1) / fs_IMU;

% % % 
% % % Definisci la lunghezza della finestra per la media mobile
% % %     window_size = 5;  % Numero di campioni
% % % 
% % %     Applica il filtro a media mobile
% % %     data_filtered = movmean(Data(4,:), window_size);
% % %     
% % % % Trova i picchi di accelerazione verticale (indicativi del contatto del tallone)
% % % [peaks, peak_locs] = findpeaks(data_filtered, 'MinPeakHeight', 0.25, 'MinPeakDistance', 0.8 * 148.15); % Imposta una distanza minima di 1.5 secondi
% % % 
% % % % Trova i minimi locali di accelerazione verticale (indicativi del sollevamento del tallone)
% % % [troughs, trough_locs] = findpeaks(-data_filtered, 'MinPeakHeight', 0.25, 'MinPeakDistance', 0.8 * 148.15); % Ricerca i minimi locali invertendo il segnale e trovando i picchi
% % % 

        signal=Data(52,:);
        [b, a] = butter(2, 10 / (fs_IMU / 2), 'low');
        filtered_signal = filtfilt(b, a, signal);
        smoothed_signal = smoothdata(filtered_signal, 'movmean', 5);
        inverted_signal = -smoothed_signal;

        [troughs, trough_locs] = findpeaks(-Data(52,:), 'MinPeakHeight', 0.25, 'MinPeakDistance', 0.8 *fs_IMU);
        [peaks, peak_locs] = findpeaks(Data(52,:), 'MinPeakHeight', 0.25, 'MinPeakDistance', 0.8 * fs_IMU); 




num_contatti_tallone = length(peaks)
num_sollevamenti_tallone = length(troughs)

% Plot dei dati e dei punti di riferimento identificati
figure;
plot(tempi_secondi,Data(52,:));
hold on;
% plot(tempi_secondi, data_filtered); hold on;
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
max_accelerazione = max(abs(Data(12,:)));
min_accelerazione = -max_accelerazione;

% Calcola i limiti dell'asse y aggiungendo un margine percentuale
margine_percentuale = 0.1; % 10% di margine
margine_superiore = max_accelerazione * (1 + margine_percentuale);
margine_inferiore = min_accelerazione * (1 - margine_percentuale);
ylim([margine_inferiore, margine_superiore]);


% % Individua gli istanti di tempo corrispondenti ai contatti del tallone e ai sollevamenti
tempi_contatti_tallone_R = tempi_secondi(peak_locs);
tempi_sollevamenti_tallone_R = tempi_secondi(trough_locs);



%% Heel's contancts extraction


clear all
close all
clc


load DATA\SUB14\SUB14_Withoutstim2.mat;

% Frequenza di campionamento
fs_IMU = 148.1481475830078; % Frequenza di campionamento in Hz

% Calcola i tempi in secondi corrispondenti ai campioni di dati
tempi_secondi = (0:length(Data(44,:))-1) / fs_IMU;

% % % 
% % % Definisci la lunghezza della finestra per la media mobile
% % %     window_size = 5;  % Numero di campioni
% % % 
% % %     Applica il filtro a media mobile
% % %     data_filtered = movmean(Data(4,:), window_size);
% % %     
% % % % Trova i picchi di accelerazione verticale (indicativi del contatto del tallone)
% % % [peaks, peak_locs] = findpeaks(data_filtered, 'MinPeakHeight', 0.25, 'MinPeakDistance', 0.8 * 148.15); % Imposta una distanza minima di 1.5 secondi
% % % 
% % % % Trova i minimi locali di accelerazione verticale (indicativi del sollevamento del tallone)
% % % [troughs, trough_locs] = findpeaks(-data_filtered, 'MinPeakHeight', 0.25, 'MinPeakDistance', 0.8 * 148.15); % Ricerca i minimi locali invertendo il segnale e trovando i picchi
% % % 

        signal=Data(44,:);
        [b, a] = butter(2, 10 / (fs_IMU / 2), 'low');
        filtered_signal = filtfilt(b, a, signal);
        smoothed_signal = smoothdata(filtered_signal, 'movmean', 5);
        inverted_signal = -smoothed_signal;

        [troughs, trough_locs] = findpeaks(-Data(44,:), 'MinPeakHeight', 0.25, 'MinPeakDistance', 0.8 *fs_IMU);
        [peaks, peak_locs] = findpeaks(Data(44,:), 'MinPeakHeight', 0.25, 'MinPeakDistance', 0.8 * fs_IMU); 




num_contatti_tallone = length(peaks)
num_sollevamenti_tallone = length(troughs)

% Plot dei dati e dei punti di riferimento identificati
figure;
plot(tempi_secondi,Data(44,:));
hold on;
% plot(tempi_secondi, data_filtered); hold on;
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
max_accelerazione = max(abs(Data(44,:)));
min_accelerazione = -max_accelerazione;

% Calcola i limiti dell'asse y aggiungendo un margine percentuale
margine_percentuale = 0.1; % 10% di margine
margine_superiore = max_accelerazione * (1 + margine_percentuale);
margine_inferiore = min_accelerazione * (1 - margine_percentuale);
ylim([margine_inferiore, margine_superiore]);


% % Individua gli istanti di tempo corrispondenti ai contatti del tallone e ai sollevamenti
tempi_contatti_tallone_L = tempi_secondi(peak_locs);
tempi_sollevamenti_tallone_L = tempi_secondi(trough_locs);



%% Plot filtered EMG data

% Frequenza di campionamento
fs = 1259.259277343750; % Frequenza di campionamento in Hz

% EMG data indices
    emg_indices = 8:8:64;

    % Initialization of the filtered signal matrix
    filtered_emg_signals = zeros(length(emg_indices), size(Data, 2));

    for i = 1:length(emg_indices)

        %Time window
        time_window_sec = 0.1;

        % Calculation of the number of samples per window(must be odd)
        window_length = round(time_window_sec * fs);
        
        if rem(window_length, 2) == 0
            window_length = window_length + 1; % Add 1 if it is even
        end

        % Calculation of the window for the moving average
        window_moving_average = ones(1, window_length) / window_length;

        % Application of the moving average using the convolution
        moving_average_signal = conv(Data(emg_indices(i),:), window_moving_average, 'same');

        % Saving the filtered signal in the filtered signals matrix
        filtered_emg_signals(i, :) = moving_average_signal;
    end

% Plot dei segnali filtrati con i nomi dei canali

for i = 1:length(emg_indices)
    figure;
    time=(0:length(filtered_emg_signals(i,:))-1)/fs;
    plot(time,filtered_emg_signals(i, :));
    xlabel('Tempo (s)');
    title(Channels(emg_indices(i),:)); % Utilizzo del nome del canale come titolo
end


%% Muscle activation gait cycle

% Definisci i nomi dei muscoli
nomi_muscoli = {'TA_R', 'VM_R', 'VL_R', 'RF_R', 'GM_R', 'ST_R', 'Gastro_R', 'SO_R'};

% Numero di muscoli
num_muscoli = numel(nomi_muscoli);

% Numero di passi
num_passi = min(length(tempi_contatti_tallone_R), length(tempi_sollevamenti_tallone_R));

% Trova la lunghezza massima di un gait cycle
lunghezza_max_ciclo = max(round((tempi_contatti_tallone_R(2:end) - tempi_contatti_tallone_R(1:end-1)) * fs));

% Creazione di un vettore normalizzato per l'asse x
x_normalized = linspace(0, 100, lunghezza_max_ciclo+1);

% Trova il valore minimo per ogni muscolo
min_values = min(filtered_emg_signals, [], 2);

% Sottrai il valore minimo da ogni segnale per normalizzarli
normalized_emg_signals = filtered_emg_signals;

% Plot dei segnali EMG normalizzati per ciascun muscolo
for m = 1:num_muscoli
    figure; % Crea una nuova figura per ciascun muscolo
    hold on; % Mantieni i grafici sullo stesso asse

    for i = 1:num_passi-16
        start_time = tempi_contatti_tallone_R(15+i);
        end_time = tempi_contatti_tallone_R(16+i);

        emg_signal = normalized_emg_signals(m, :);

        % Converti i tempi in secondi in indici dei campioni
        start_index = round(start_time * fs); % fs è la frequenza di campionamento del segnale EMG
        end_index = round(end_time * fs);

        % Interpola il segnale EMG per avere la stessa lunghezza per ogni gait cycle
        emg_interpolato = interp1(start_index:end_index, emg_signal(start_index:end_index), ...
            linspace(start_index, end_index, lunghezza_max_ciclo+1));

% %         % Plot del segnale EMG del muscolo corrente per il ciclo corrente
%         plot(x_normalized, emg_interpolato, 'LineWidth', 1.5);

        mean_values_R(i,m)=mean(emg_interpolato);

        % Memorizza il segnale interpolato nella matrice dei cicli EMG
        emg_cicli(:,i) = emg_interpolato';
    end
    
    % Calcola la media dei segnali EMG per il muscolo corrente su tutti i cicli di camminata
    emg_media = mean(emg_cicli, 2);

     % Calcola la deviazione standard dei segnali EMG per il muscolo corrente su tutti i cicli di camminata
    emg_std = std(emg_cicli, 0, 2); % Calcola la deviazione standard lungo le colonne (dimensione 2)


    % Plot della media dei segnali EMG sopra i singoli grafici
    plot(x_normalized, emg_media, 'k', 'LineWidth', 5); % 'k' indica il colore nero
% 
%         % Plot dell'intervallo di confidenza della deviazione standard
%     plot(x_normalized, emg_media + 2*emg_std, '--r', 'LineWidth', 1); % +2 STD
%     plot(x_normalized, emg_media - 2*emg_std, '--r', 'LineWidth', 1); % -2 STD

    title([ nomi_muscoli{m}]);
    xlabel('Gait Cycle (%)');
    ylabel('Activation percentage');
    hold off; % Termina il mantenimento dei grafici sullo stesso asse
end


%% detection
   % Conversione dell'intervallo di tempo in indici di campionamento
%     start_index = round(tempi_contatti_tallone_R(30) * fs_IMU);
%     end_index = round(tempi_contatti_tallone_R(31) * fs_IMU);



start_index=peak_locs(1);
end_index=peak_locs(2);
    % Estrazione del segnale nell'intervallo specificato
    feet_acc = (Data(12,start_index:end_index)/9.81);

    % Calcolo di dt
    dt = 1 / fs_IMU;

    % Integrazione dei dati del giroscopio filtrati per ottenere gli angoli
    feet_vel = cumtrapz((0:length(feet_acc)-1)*dt, feet_acc); 
    feet_pos= cumtrapz((0:length(feet_vel)-1)*dt, feet_vel); 

    figure;
    subplot(2,1,1);
    plot(feet_vel)
    subplot(2,1,2);
    plot(feet_pos)