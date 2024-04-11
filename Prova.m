clear all;
close all;
clc;


load Gait\50HZ\50HZ_RF.mat;

emg_data1=emg_data(50000:100000);

x=abs(emg_data1);
% Calcola il filtro RMS
finestra = 840; % Dimensione della finestra per il calcolo RMS
segnale_rms = sqrt(movmean(x.^2, finestra)); % Calcolo del filtro RMS

% Parametri del filtro passa banda Butterworth
lowcut = 20; % Frequenza di taglio bassa (Hz)
highcut = 150; % Frequenza di taglio alta (Hz)
order_butter = 6; % Ordine del filtro Butterworth

% Parametri del filtro a doppio T-trap
f0 = 100; % Frequenza di rete elettrica (Hz)
order_ttrap = 2; % Ordine del filtro a T

% Frequenza di campionamento
fs = 2500; % Hz

% Applicazione del filtro passa banda Butterworth
[b_butter, a_butter] = butter(order_butter, [lowcut, highcut]/(fs/2), 'stop');
x_filtered_butter = filtfilt(b_butter, a_butter, emg_data1); % x è il segnale di ingresso

% Frequenze normalizzate per i filtri a T
fn1 = f0 / (fs / 2); % Frequenza normalizzata per il primo filtro a T
fn2 = f0 / (fs / 2); % Frequenza normalizzata per il secondo filtro a T

% Progettazione dei filtri a T
[b1, a1] = butter(order_ttrap, [fn1-0.02 fn1+0.02], 'stop'); % Primo filtro a T
[b2, a2] = butter(order_ttrap, [fn2-0.02 fn2+0.02], 'stop'); % Secondo filtro a T

% Combinazione dei filtri a T per formare il filtro a doppio T-trap
b_double_T_trap = conv(b1, b2);
a_double_T_trap = conv(a1, a2);

% Applicazione del filtro a doppio T-trap
x_filtered_combined = filtfilt(b_double_T_trap, a_double_T_trap, x_filtered_butter);

% Applicazione del denoising wavelet
livello = 4; % Livello di decomposizione wavelet
segnale_filtrato = wdenoise(x_filtered_combined, livello);

Rectified_data = filtfilt(ones(1,840)/840,1,abs(segnale_filtrato));

t=1:1:length(emg_data1);

% Plot del segnale originale e del segnale filtrato
figure;
plot(t, emg_data1, 'b-', 'LineWidth', 1.5); % Plot del segnale originale
hold on;
plot(t, Rectified_data, 'r-', 'LineWidth', 1.5); % Plot del segnale filtrato
xlabel('Tempo [s]');
ylabel('Ampiezza');
title('Filtro composto (Passa banda Butterworth + Doppio T-trap)');
legend('Segnale originale', 'Segnale filtrato');
grid on;

figure;
subplot(3,1,1);
plot(emg_data1); hold on;
subplot(3,1,2)
plot(segnale_rms); hold on;
subplot(3,1,3);
plot(Rectified_data); hold on;



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

