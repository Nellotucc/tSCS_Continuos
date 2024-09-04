% Definisci il percorso completo del file .fig
filePath = 'PLOT\SUB1\ST_DIST_L.fig';  % Modifica questo percorso con il tuo

% Carica il file .fig
fig = openfig(filePath, 'invisible');

% Trova tutti gli assi (subplot) nella figura
allAxes = findall(fig, 'Type', 'Axes');  % Trova tutti gli assi nella figura

% Prepara una struttura per memorizzare i dati
dataStruct = struct();

% Estrai i dati da ciascun subplot
for k = 1:length(allAxes)
    % Trova tutti gli oggetti di tipo 'line' all'interno dell'asse corrente
    h = findobj(allAxes(k), 'Type', 'line');
    
    % Inizializza celle per memorizzare i dati di questo subplot
    subplotData = struct();
    
    % Estrai i dati dai grafici
    for i = 1:length(h)
        xData = get(h(i), 'XData');
        yData = get(h(i), 'YData');
        
        % Memorizza i dati nella struttura
        subplotData(i).xData = xData;
        subplotData(i).yData = yData;
    end
    
    % Aggiungi i dati di questo subplot alla struttura principale
    dataStruct.(sprintf('subplot_%d', k)) = subplotData;
end

% Salva la struttura di dati in un file .mat
save('all_subplots_data.mat', 'dataStruct');

% Chiudi la figura
close(fig);
