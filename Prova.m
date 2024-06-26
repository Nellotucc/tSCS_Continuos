
%% Load segnale 
clear all;
close all;
clc;


load Gait\50HZ\50HZ_SOL.mat;

emg_data1=emg_data(50000:100000);

x=abs(emg_data1);


%% rectification and zero-lag with 3hz cutoff
cutoff_frequency = 3; % Frequenza di taglio desiderata in Hz
fs = 2500; % Frequenza di campionamento in Hz
num_samples = round(fs / cutoff_frequency); % Calcolo del numero di campioni per la media mobile

% Coefficienti del filtro
b = ones(1, num_samples) / num_samples;
a = 1;

% Applicazione del filtro utilizzando filtfilt
Rectified_data = filtfilt(b, a, x);


%% Low pass butterworth 3 hz


ordine = 5; % Ordine del filtro
frequenza_di_taglio = 3; % Frequenza di taglio in Hz
frequenza_di_campionamento = 2500; % Frequenza di campionamento in Hz

% Calcola le frequenze normalizzate
frequenza_normalizzata = frequenza_di_taglio / (frequenza_di_campionamento / 2);

% Calcola i coefficienti del filtro di Butterworth
[b, a] = butter(ordine, frequenza_normalizzata, 'low');

% Applica il filtro utilizzando filtfilt
butterworth = filtfilt(b, a, x);



%% Band pass 40-1000Hz + notch 60Hz 
% Parametri del filtro passa banda
frequenza_passa_banda = [40 1000]; % Frequenze di taglio della banda passante in Hz

% Calcola le frequenze normalizzate
frequenze_normalizzate_passa_banda = frequenza_passa_banda / (frequenza_di_campionamento / 2);

% Calcola i coefficienti del filtro passa banda di Butterworth
[b_passa_banda, a_passa_banda] = butter(5, frequenze_normalizzate_passa_banda, 'bandpass');

segnale_passa_banda = filtfilt(b_passa_banda, a_passa_banda, emg_data1);

% Parametri del filtro notch
frequenza_notch = 60; % Frequenza del filtro notch in Hz

% Calcola la frequenza normalizzata per il filtro notch
frequenza_normalizzata_notch = frequenza_notch / (frequenza_di_campionamento / 2);

% Calcola i coefficienti del filtro notch di Butterworth
[b_notch, a_notch] = butter(5, [frequenza_normalizzata_notch-0.02, frequenza_normalizzata_notch+0.02], 'stop');

% Applica il filtro passa banda seguito dal filtro notch utilizzando filtfilt
bandandnotch = filtfilt(b_notch, a_notch, filtfilt(b_passa_banda, a_passa_banda, segnale_passa_banda));

abs_bandnotch=abs(bandandnotch);

% Applicazione del filtro utilizzando filtfilt
xx = filtfilt(ones(1, num_samples) / num_samples, 1, bandandnotch);

%%  Plot
t=1:1:length(emg_data1);
figure;
plot(t, emg_data1, 'b-', 'LineWidth', 1.5); % Plot del segnale originale
hold on;
plot(t, Rectified_data, 'r-', 'LineWidth', 1.5); % Plot del segnale filtrato
xlabel('Tempo [s]');
ylabel('Ampiezza');
title('Filtro composto (Passa banda Butterworth + Doppio T-trap)');
legend('Segnale originale', 'Segnale filtrato');
grid on;



%% Plot
figure;
subplot(5,1,1);
plot(emg_data1); hold on;
subplot(5,1,2)
plot(x); hold on;
subplot(5,1,3);
plot(Rectified_data); hold on;
subplot(5,1,4);
plot(butterworth); hold on;
subplot(5,1,5);
plot(xx); hold on;


%% Confronto
% Calcola la frequenza di campionamento
frequenza_campionamento = 2500;

% Definisci il vettore tempo
tempo = (0:length(xx)-1) / frequenza_campionamento;

% Normalizza i segnali
butterworth_normalized = butterworth / max(abs(butterworth));
xx_normalized = xx / max(abs(xx));
Rectified_data_normalized = Rectified_data / max(abs(Rectified_data));

% Definisci l'offset per ciascun segnale
offset_butterworth = 2;
offset_xx = 1;
offset_rectified_data = 0;

% Plotta i segnali shiftati verticalmente
figure;
plot(tempo, butterworth_normalized + offset_butterworth, 'r', 'LineWidth', 2);
hold on;
plot(tempo, xx_normalized + offset_xx, 'g', 'LineWidth', 2);
plot(tempo, Rectified_data_normalized + offset_rectified_data, 'b', 'LineWidth', 2);

% Aggiungi etichette e legenda
xlabel('Tempo (s)');
ylabel('Valore Normalizzato del Segnale');
title('Confronto dei Segnali Normalizzati e Shiftati Verticalmente');
legend('rect+lowpass', 'bandpass+notch+rect', 'Rect+zerolag');



% Rectified_data = filtfilt(ones(1,20)/20,1,abs(emg_data1));
% y = fft(Rectified_data);
% y(1) = [];
% figure;
% plot(y,'ro')
% xlabel('real(y)')
% ylabel('imag(y)')
% title('Fourier Coefficients')
% 
% 
% 
% n = length(Rectified_data);
% power = abs(Rectified_data(1:floor(n/2))).^2; % power of first half of transform data
% maxfreq = 1/2;                   % maximum frequency
% freq = (1:n/2)/(n/2)*maxfreq;
% figure; % equally spaced frequency grid
% plot(freq,power)
% 

