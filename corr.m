% Definisci il percorso della cartella principale contenente le sottocartelle
mainFolderPath = 'PLOT';  % Modifica questo percorso con il tuo

% Ottieni una lista di tutte le sottocartelle nella cartella principale
subFolders = dir(mainFolderPath);
subFolders = subFolders([subFolders.isdir] & ~ismember({subFolders.name}, {'.', '..'}));

% Prepara una cella per memorizzare i dati da tutti i file e sottocartelle
allData = {};

% Cicla attraverso ciascuna sottocartella
for folderIdx = 1:length(subFolders)
    % Ottieni il percorso completo della sottocartella
    folderPath = fullfile(mainFolderPath, subFolders(folderIdx).name);
    
    % Ottieni una lista di tutti i file .fig nella sottocartella
    fileList = dir(fullfile(folderPath, '*.fig'));
    
    % Cicla attraverso ciascun file .fig nella sottocartella
    for fileIdx = 1:length(fileList)
        % Ottieni il percorso completo del file .fig
        filePath = fullfile(folderPath, fileList(fileIdx).name);
        
        % Carica il file .fig
        fig = openfig(filePath, 'invisible');
        
        % Trova tutti gli assi (subplot) nella figura
        allAxes = findall(fig, 'Type', 'Axes');
        
        % Cicla attraverso ciascun subplot e memorizza i dati
        for k = 1:length(allAxes)
            % Trova tutti gli oggetti di tipo 'line' all'interno dell'asse corrente
            h = findobj(allAxes(k), 'Type', 'line');
            
            % Estrai i dati dai grafici
            for i = 1:length(h)
                xData = get(h(i), 'XData');
                yData = get(h(i), 'YData');
                
                % Salva i dati estratti in una cella
                % Ogni riga della cella rappresenta un grafico
                % Ogni colonna rappresenta: [nome_sottocartella, nome_file, subplotIdx, xData, yData]
                allData{end+1, 1} = subFolders(folderIdx).name;  % Nome della sottocartella
                allData{end, 2} = fileList(fileIdx).name;  % Nome del file
                allData{end, 3} = k;  % Indice del subplot
                allData{end, 4} = xData;  % Dati x
                allData{end, 5} = yData;  % Dati y
            end
        end
        
        % Chiudi la figura
        close(fig);
    end
end

% Salva i dati in un file .mat
save('all_subplots_data.mat', 'allData');

