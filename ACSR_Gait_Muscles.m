
%% filtraggio
clear all;
close all;
clc;

% Percorso della directory principale
mainPath = 'DATA\Gait\Muscles';

% Ottieni la lista delle sottodirectory nella directory principale
subdirectories = dir(mainPath);
subdirectories = subdirectories([subdirectories.isdir]); % Filtra solo le directory
subdirectories = subdirectories(~ismember({subdirectories.name}, {'.', '..'})); % Rimuovi le cartelle speciali . e ..

% Ciclo attraverso le sottodirectory
for dirIdx = 1:numel(subdirectories)
    currentDir = fullfile(mainPath, subdirectories(dirIdx).name);
    
    % Ottieni la lista dei file .mat nella sottodirectory corrente
    files = dir(fullfile(currentDir, '*.mat'));

figure;
for i = 1:numel(files)
    File_name = fullfile(currentDir, files(i).name);% Obtain current file name
        load(File_name);  

 %Filter 
 emg_data1=(emg_data)';
 stop=uint32(length(emg_data1)/5);
 emg_for_training1=emg_data1(1,1:37500);
 emg_for_training2=emg_data1(1,45000:100000);
 ACSR_window=800;

emg_filtered1=ACSR_filter(emg_for_training1,emg_data1,ACSR_window);
emg_filtered2=ACSR_filter(emg_for_training2,emg_data1,ACSR_window);

    [~, name, ~] = fileparts(File_name);  %extract name file
     titolo = strrep(name, '_', ' ');  % CHANGE THE UNDERSCORE WITH SPACE!


  time=[1:1:length(emg_data1)];

   subplot(2,numel(files),i);
    plot(time,emg_data1,'b');hold on;
     xlabel('Time [s]');ylabel('Amplitude [mV]');
    title(titolo,'fontsize',12,'fontweight','bold');

subplot(2,numel(files),i+numel(files));
    plot(time,emg_filtered2,'r');
    xlabel('Time [s]');ylabel('Amplitude [mV]');
    title('Artifact filtered','fontsize',12,'fontweight','bold');
  end


end
