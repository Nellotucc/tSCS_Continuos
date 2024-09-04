% Estrazione ciclo di camminata da dati IMU

clear all
close all
clc


load DATA\03_06_AHMED\AHMED_Withoutstim_1.mat;

%% Cambio delle prime 8 righe 
% Supponiamo che la tua matrice si chiami Data
% Estrai le dimensioni della matrice
[nRows, nCols] = size(Data);

% Controlla se ci sono almeno 8 righe
if nRows >= 8
    % Estrai le prime 8 righe
    firstEightRows = Data(1:8, :);
    
    % Estrai le restanti righe
    remainingRows = Data(9:end, :);
    
    % Concatenale per ottenere la nuova matrice
    newData = [remainingRows; firstEightRows];
else
    error('La matrice non ha almeno 8 righe.');
end

% Assegna la nuova matrice a Data
Data = newData;

%% Heel's contancts extraction

% Frequenza di campionamento
fs_IMU = 148.1481475830078; % Frequenza di campionamento in Hz

% Calcola i tempi in secondi corrispondenti ai campioni di dati
tempi_secondi = (0:length(Data(4,:))-1) / fs_IMU;

% Definisci la lunghezza della finestra per la media mobile
    window_size = 12;  % Numero di campioni

    % Applica il filtro a media mobile
    data_filtered = movmean(Data(4,:), window_size);
    

% Trova i picchi di accelerazione verticale (indicativi del contatto del tallone)
[peaks, peak_locs] = findpeaks(Data(4,:), 'MinPeakHeight', 0.25, 'MinPeakDistance', 0.8 * 148.15); % Imposta una distanza minima di 1.5 secondi

% Trova i minimi locali di accelerazione verticale (indicativi del sollevamento del tallone)
[troughs, trough_locs] = findpeaks(-data_filtered, 'MinPeakHeight', 0.25, 'MinPeakDistance', 0.8 * 148.15); % Ricerca i minimi locali invertendo il segnale e trovando i picchi

num_contatti_tallone = length(peaks)
num_sollevamenti_tallone = length(troughs)

% Plot dei dati e dei punti di riferimento identificati
figure;
plot(tempi_secondi, Data(4,:)); hold on;
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
max_accelerazione = max(abs(Data(68,:)));
min_accelerazione = -max_accelerazione;

% Calcola i limiti dell'asse y aggiungendo un margine percentuale
margine_percentuale = 0.1; % 10% di margine
margine_superiore = max_accelerazione * (1 + margine_percentuale);
margine_inferiore = min_accelerazione * (1 - margine_percentuale);
ylim([margine_inferiore, margine_superiore]);


% % Individua gli istanti di tempo corrispondenti ai contatti del tallone e ai sollevamenti
tempi_contatti_tallone_R = tempi_secondi(peak_locs);
tempi_sollevamenti_tallone_R = tempi_secondi(trough_locs);



%% Plot filtered EMG data


% Indici delle righe EMG
emg_indices = 1:8:size(Data,1);

% Frequenza di campionamento
fs = 1259.259277343750; % Frequenza di campionamento in Hz

% Parametri del filtro passa banda Butterworth
fc_low = 10; % Frequenza di taglio inferiore in Hz
fc_high = 500; % Frequenza di taglio superiore in Hz
order_bandpass = 4; % Ordine del filtro passa banda

% Parametri del filtro passa basso Butterworth (per la media mobile)
fc_lp = 3; % Frequenza di taglio del filtro passa basso in Hz
order_lowpass = 4; % Ordine del filtro passa basso

% Creazione del filtro passa banda
[b_bandpass, a_bandpass] = butter(order_bandpass, [fc_low/(fs/2), fc_high/(fs/2)], 'bandpass');

% Creazione del filtro passa basso (per la media mobile)
[b_lowpass, a_lowpass] = butter(order_lowpass, fc_lp/(fs/2), 'low');

for i = 1:length(emg_indices)

    % Applicazione del filtro passa banda
    emg_filtered(i,:) = filtfilt(b_bandpass, a_bandpass, Data(emg_indices(i), :));
    
    % Raddrizzamento dell'onda
    emg_rectified(i,:) = abs(emg_filtered(i,:));
    % Converti la finestra temporale da millisecondi a secondi
    finestra_temporale_sec = 0.1;

% Calcola il numero di campioni nella finestra (deve essere dispari)
lunghezza_finestra = round(finestra_temporale_sec * fs);
% floor(lunghezza_finestra / 2);

if rem(lunghezza_finestra, 2) == 0
    lunghezza_finestra = lunghezza_finestra + 1; % Aggiungi 1 se è pari
end

% Calcola la finestra di media mobile
finestra_media_mobile = ones(1, lunghezza_finestra) / lunghezza_finestra;

% Applica la media mobile al segnale utilizzando la convoluzione
segnale_media_mobile = conv(emg_rectified(i,:), finestra_media_mobile, 'same');
    
%     % Applicazione del filtro passa basso per la media mobile
%     emg_smoothed = filtfilt(b_lowpass, a_lowpass, emg_rectified(i,:));


    % Salvataggio del segnale filtrato nella matrice dei segnali filtrati
    filtered_emg_signals(i, :) = segnale_media_mobile;
end


% Plot dei segnali filtrati con i nomi dei canali

for i = 1:length(emg_indices)
    figure;
    time=(0:length(filtered_emg_signals(i,:))-1)/fs;
    plot(time,filtered_emg_signals(i, :),'k');
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
normalized_emg_signals = filtered_emg_signals - min_values;

% Plot dei segnali EMG normalizzati per ciascun muscolo
for m = 1:num_muscoli
    figure; % Crea una nuova figura per ciascun muscolo
    hold on; % Mantieni i grafici sullo stesso asse

    for i = 1:num_passi-1
        start_time = tempi_contatti_tallone_R(i);
        end_time = tempi_contatti_tallone_R(i+1);

        emg_signal = normalized_emg_signals(m, :);

        % Converti i tempi in secondi in indici dei campioni
        start_index = round(start_time * fs); % fs è la frequenza di campionamento del segnale EMG
        end_index = round(end_time * fs);

        % Interpola il segnale EMG per avere la stessa lunghezza per ogni gait cycle
        emg_interpolato = interp1(start_index:end_index, emg_signal(start_index:end_index), ...
            linspace(start_index, end_index, lunghezza_max_ciclo+1));

%         % Plot del segnale EMG del muscolo corrente per il ciclo corrente
        plot(x_normalized, emg_interpolato, 'LineWidth', 1.5);

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

    title(['Segnale EMG - ', nomi_muscoli{m}]);
    xlabel('Gait Cycle (%)');
    ylabel('Amplitude');
    hold off; % Termina il mantenimento dei grafici sullo stesso asse
end
