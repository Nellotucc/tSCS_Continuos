clear all
close all
clc


<<<<<<< HEAD
Path = 'C:\Users\39320\Dropbox (Politecnico Di Torino Studenti)\PC\Desktop\EPFL\Thesis code\tSCS_TNE\MATLAB\DATA\DOME\30Hz\L1-L2';
=======
Path = 'C:\Users\local_B216353\Documents\tSCS\tSCS_TNE\MATLAB\DATA\DOME\30Hz\L1-L2';
>>>>>>> 53f70bf10cb2d95fdb1b6068a525d8679b66d1e0


files = dir(fullfile(Path, '*.mat'));  % File list

for i = 1:numel(files)
    File_name = fullfile(Path, files(i).name); % Obtain current file name
    load(File_name); 


    %% Filter 
 stop=uint32(length(emg_data)*5/6);
 emg_for_training1=emg_data(1,1:4000);
 emg_for_training2=emg_data(1,1:stop);
 ACSR_window=200;

emg_filtered1=ACSR_filter(emg_for_training1,emg_data,ACSR_window);
emg_filtered2=ACSR_filter(emg_for_training2,emg_data,ACSR_window);

    [~, name, ~] = fileparts(File_name);  %extract name file
     titolo = strrep(name, '_', ' ');  % CHANGE THE UNDERSCORE WITH SPACE!

   figure;
  time=[1:1:length(emg_data)];
   subplot(2,2,1);
    plot(time,emg_data,'b');hold on;
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
pos(3) = 0.5;  % width
pos(4) = 0.35; % height 
set(subplot(2,2,1), 'Position', pos);

end