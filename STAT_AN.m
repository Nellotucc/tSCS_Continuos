clear all
close all
clc



% Definisci la cartella contenente i file .mat
folder = 'MATRIX\HIP_STANCE';
files = dir(fullfile(folder, '*.mat'));

% Inizializza una variabile per la concatenazione
allData = [];

% Cicla attraverso ciascun file
for k = 1:length(files)
    % Ottieni il percorso completo del file
    filePath = fullfile(folder, files(k).name);
    
    % Carica il file .mat
    data = load(filePath);
    
    % Supponiamo che il nome della variabile nel file sia 'matrix'
    % (Sostituisci 'matrix' con il nome della variabile se diverso)
    fieldName = fieldnames(data);
    if length(fieldName) ~= 1
        error('File %s contiene più di una variabile. Verifica il contenuto del file.', files(k).name);
    end
    matrix = data.(fieldName{1});
    
    % Concatenare orizzontalmente la matrice
    allData = [allData, matrix];
end

% Verifica la dimensione della matrice concatenata
disp(size(allData));

% Salva la matrice concatenata in un nuovo file .mat se necessario
save('ConcatenatedData.mat', 'allData');

% Supponiamo che allData sia la tua matrice
[numRows, numCols] = size(allData);

% Inizializza una matrice con NaN (o un altro valore di riempimento)
maxLength = numCols; % Puoi scegliere un valore che rappresenta una lunghezza massima
filteredData = NaN(numRows, maxLength);

% Cicla attraverso ciascuna riga
for i = 1:numRows
    % Estrai la riga i-esima
    row = allData(i, :);
    
    % Rimuovi gli zeri dalla riga
    filteredRow = row(row ~= 0);
    
    % Trova la lunghezza della riga filtrata
    rowLength = length(filteredRow);
    
    % Inserisci la riga filtrata nella matrice
    filteredData(i, 1:rowLength) = filteredRow;
end

% Verifica la dimensione della matrice finale
disp(size(filteredData));

% Salva la matrice senza zeri in un nuovo file .mat
save('FilteredData.mat', 'filteredData');



frequencies = [100,15,30,50,5,65,75,0];
% Sorting frequencies and relative percentage
[sorted_frequencies, sort_idx] = sort(frequencies);

filteredData_PROX1=filteredData(8:15,:);
filteredData_PROX=filteredData_PROX1(sort_idx,:);
% Supponiamo che filteredData sia la tua matrice 16x1551
% Seleziona le prime 7 righe
first7Rows = filteredData(1:7, :);

% Seleziona la 15ª riga
row15 = filteredData(15, :);

% Combina le prime 7 righe con la 15ª riga
% Per assicurarti che le righe abbiano la stessa dimensione, concatena verticalmente
filteredData_DIST1 = [first7Rows; row15];
filteredData_DIST=filteredData_DIST1(sort_idx,:);


%% stat DIST

[numRighe, numColonne] = size(filteredData_DIST);
stop=numColonne-100;

group1 = filteredData_DIST(1, 5:stop);
group2 = filteredData_DIST(2, 5:stop);
group3 = filteredData_DIST(3, 5:stop);
group4 = filteredData_DIST(4, 5:stop);
group5 = filteredData_DIST(5, 5:stop);
group6 = filteredData_DIST(6, 5:stop);
group7 = filteredData_DIST(7, 5:stop);
group8 = filteredData_DIST(8, 5:stop);


% Creation of a unic group
data = [group1, group2, group3, group4, group5, group6, group7, group8];

% Frequencies associatated to each group
frequencies = {'0Hz', '5Hz', '15Hz', '30Hz', '50Hz','65Hz', '75Hz', '100Hz'};

%Global group
group = [repmat(frequencies(1), 1, length(group1)), ...
         repmat(frequencies(2), 1, length(group2)), ...
         repmat(frequencies(3), 1, length(group3)), ...
         repmat(frequencies(4), 1, length(group4)), ...
         repmat(frequencies(5), 1, length(group5)), ...
         repmat(frequencies(6), 1, length(group6)), ...
         repmat(frequencies(7), 1, length(group7)), ...
         repmat(frequencies(8), 1, length(group8))];
      

%Unidirectional ANOVA
[p, tbl, stats] = anova1(data, group);
    %titleText = inputname(1);
    %title(titleText);

% ANOVA Results
disp('Tabella ANOVA:');
disp(tbl);

figure;
% Test Newman-Keuls (Tukey-Kramer in MATLAB)
results = multcompare(stats, 'CType', 'tukey-kramer');
disp('Risultati del test post-hoc:');
disp(results);


%% stat prox

[numRighe, numColonne] = size(filteredData_PROX);
stop=numColonne-100;

group1 = filteredData_PROX(1, 5:stop);
group2 = filteredData_PROX(2, 5:stop);
group3 = filteredData_PROX(3, 5:stop);
group4 = filteredData_PROX(4, 5:stop);
group5 = filteredData_PROX(5, 5:stop);
group6 = filteredData_PROX(6, 5:stop);
group7 = filteredData_PROX(7, 5:stop);
group8 = filteredData_PROX(8, 5:stop);


% Creation of a unic group
data = [group1, group2, group3, group4, group5, group6, group7, group8];

% Frequencies associatated to each group
frequencies = {'0Hz', '5Hz', '15Hz', '30Hz', '50Hz','65Hz', '75Hz', '100Hz'};

%Global group
group = [repmat(frequencies(1), 1, length(group1)), ...
         repmat(frequencies(2), 1, length(group2)), ...
         repmat(frequencies(3), 1, length(group3)), ...
         repmat(frequencies(4), 1, length(group4)), ...
         repmat(frequencies(5), 1, length(group5)), ...
         repmat(frequencies(6), 1, length(group6)), ...
         repmat(frequencies(7), 1, length(group7)), ...
         repmat(frequencies(8), 1, length(group8))];
      

%Unidirectional ANOVA
[p, tbl, stats] = anova1(data, group);
    %titleText = inputname(1);
    %title(titleText);

% ANOVA Results
disp('Tabella ANOVA:');
disp(tbl);

figure;
% Test Newman-Keuls (Tukey-Kramer in MATLAB)
results = multcompare(stats, 'CType', 'tukey-kramer');
disp('Risultati del test post-hoc:');
disp(results);