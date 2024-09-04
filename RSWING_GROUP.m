
clear all
close all
clc


% Definisci il percorso della cartella
data_dir = 'C:\Users\39320\Desktop\EPFL\Thesis code\tSCS_Continuos\NEW_PROTOCOL\MATRIX\EMG\RSWING'; 

% Ottieni tutti i file .mat nella directory
mat_files = dir(fullfile(data_dir, '*.mat'));

% Mostra i file trovati
disp('File trovati:');
disp({mat_files.name});  % Mostra i nomi dei file

% Verifica il numero di file
numFiles = length(mat_files);
if numFiles == 0
    error('Nessun file .mat trovato nella cartella: %s', data_dir);
end
% Inizializza una matrice di zeri con le dimensioni delle matrici da sommare (7x16)
sumMatrix = zeros(7, 16);

% Itera sui file e carica le matrici
for i = 1:numFiles
    % Crea il nome del file
    fileName = mat_files(i).name;
    
    % Crea il percorso completo
    fullPath = fullfile(data_dir, fileName);
    
    % Verifica se il file esiste
    if exist(fullPath, 'file')
        % Carica la matrice dal file
        tempStruct = load(fullPath);  % Carica la matrice in una struttura temporanea
        
        % Estrai il campo della variabile
        fieldNames = fieldnames(tempStruct);
        tempMatrix = tempStruct.(fieldNames{1});  % Usa il primo nome di campo
        
        % Sostituisci NaN con zeri nella matrice temporanea
        tempMatrix(isnan(tempMatrix)) = 0;
        
        % Loop su ogni riga della matrice
        for j = 1:size(tempMatrix, 1)
            % Controlla se almeno un elemento della riga è superiore a 100
            if any(tempMatrix(j, :) > 100)
                % Trova il massimo valore nella riga
                max_value = max(tempMatrix(j, :));
                
                % Normalizza la riga dividendo ogni elemento per il massimo
                tempMatrix(j, :) = (tempMatrix(j, :) / max_value) * 100;
            end
        end
        
        % Somma la matrice normalizzata alla matrice somma
        sumMatrix = sumMatrix + tempMatrix;
    else
        error(['File non trovato: ', fullPath]);
    end
end

disp('Somma delle matrici completata.');



sumMatrixTOT=(sumMatrix)/12;

muscles_R = {'GM_R', 'GASTRO_R','TA_R', 'ST_R', 'VM_R', 'VL_R', 'RF_R'};

% Frequencies of the corrisponding value
frequencies = [100,15,30,50,5,65,75,0];

% Delete the last column
 max_values_RswingPROX = sumMatrixTOT(:, [1:7,15]);
 max_values_RswingDIST = sumMatrixTOT(:, 8:15);

 % Sort the frequencies and relative peak values
 [sorted_frequencies, sort_idx] = sort(frequencies);
 
 max_values_Rswing_sortedPROX = max_values_RswingPROX(:, sort_idx);
 max_values_Rswing_sortedDIST = max_values_RswingDIST(:, sort_idx);


 % Savitzky-Golay Filter parameters
sgolay_order = 1; % Polynomial order
sgolay_frame = 3; % Window length (bigger, bigger smooth action

for i = 1:7
    subplot(2, 4, i); 
    
    % Apply Savitzky-Golay filter
    smoothed_values_stance = sgolayfilt(max_values_Rswing_sortedPROX(i, :), sgolay_order, sgolay_frame);
    %semilogx(sorted_frequencies,smoothed_values_stance, '-o', 'LineWidth', 2); % Plot logaritmico
    plot(sorted_frequencies,smoothed_values_stance, 'LineWidth', 2); % Plot logaritmico
    hold on;
    title('Swing phase PROX', muscles_R{i});
    xlabel('Frequency (Hz)');
    ylabel('Peak Value');
    %grid on;
    hold off;
end

figure;
for i = 1:7
    subplot(2, 4, i); 
    
    % Apply Savitzky-Golay filter
    smoothed_values_stance = sgolayfilt(max_values_Rswing_sortedDIST(i, :), sgolay_order, sgolay_frame);
    %semilogx(sorted_frequencies,smoothed_values_stance, '-o', 'LineWidth', 2); % Plot logaritmico
    plot(sorted_frequencies,smoothed_values_stance, 'LineWidth', 2); % Plot logaritmic
    hold on;
%     xline(100, '--r', 'LineWidth', 1.5); % Linea verticale tratteggiata rossa a 100 Hz
    title('Swing phase DIST', muscles_R{i});
    xlabel('Frequency (Hz)');
    ylabel('Peak Value');
    grid on;
    hold off;
end

