
%% conversione
clear all;
close all;
clc;

% Specifica la directory di salvataggio
directory = 'DATA\Gait\30HZ';

% Ottieni la lista di tutti i file .txt nella directory
filelist = dir(fullfile(directory, '*.txt'));

% Ciclo for per attraversare tutti i file .txt nella cartella
for i = 1:length(filelist)
    % Leggi i dati dal file .txt corrente
    data = importdata(fullfile(directory, filelist(i).name));
    
    % Estrai i dati della corrente (Channel 1)
    emg_data = data.data;
    
    % Crea il nome del file .mat da salvare
    [~, filename, ~] = fileparts(filelist(i).name);
    mat_filename = fullfile(directory, [filename '.mat']);
    
    % Salva i dati in un file .mat nella stessa cartella
    save(mat_filename, 'emg_data');
end


pause(5);

%% filtraggio
clear all;
close all;
clc;


Path = 'DATA\Gait\50HZ';


files = dir(fullfile(Path, '*.mat'));  % File list

for i = 1:numel(files)
    File_name = fullfile(Path, files(i).name);% Obtain current file name
        load(File_name);  

 %Filter 
 emg_data1=(emg_data)';
 stop=uint32(length(emg_data1)/5);
 emg_for_training1=emg_data1(1,1:37500);
 emg_for_training2=emg_data1(1,37500:100000);
 ACSR_window=800;

emg_filtered1=ACSR_filter(emg_for_training1,emg_data1,ACSR_window);
emg_filtered2=ACSR_filter(emg_for_training2,emg_data1,ACSR_window);

    [~, name, ~] = fileparts(File_name);  %extract name file
     titolo = strrep(name, '_', ' ');  % CHANGE THE UNDERSCORE WITH SPACE!

   figure;
  time=[1:1:length(emg_data1)];
   subplot(2,2,1);
    plot(time,emg_data1,'b');hold on;
     xlabel('Time [s]');ylabel('Amplitude [mV]');
    title(titolo,'fontsize',12,'fontweight','bold');
subplot(2,2,3);
    plot(time,emg_filtered1,'r');
    xlabel('Time [s]');ylabel('Amplitude [mV]');
    title('Noise-filtered','fontsize',12,'fontweight','bold');
subplot(2,2,4);
    plot(time,emg_filtered2,'r');
    xlabel('Time [s]');ylabel('Amplitude [mV]');
    title('Noise&Artifact filtered','fontsize',12,'fontweight','bold');

pos = get(subplot(2,2,1), 'Position');
pos(1) = 0.25; % move right
pos(2) = 0.6;  % Move up
pos(3) = 0.5;  % change width
pos(4) = 0.35; % change height 
set(subplot(2,2,1), 'Position', pos);

end





