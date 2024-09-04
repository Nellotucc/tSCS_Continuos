% Inizializzazione
figurePaths = {'PLOT/SUB1/SW_DIST_R.fig', ...
               'PLOT/SUB2/SW_DIST_R.fig', ...
               'PLOT/SUB3/SW_DIST_R.fig', ...
               'PLOT/SUB4/SW_DIST_R.fig', ...
               'PLOT/SUB5/SW_DIST_R.fig', ...
               'PLOT/SUB6/SW_DIST_R.fig', ...
               'PLOT/SUB7/SW_DIST_R.fig', ...
               'PLOT/SUB8/SW_DIST_R.fig', ...
               'PLOT/SUB9/SW_DIST_R.fig', ...
               'PLOT/SUB10/SW_DIST_R.fig', ...
               'PLOT/SUB11/SW_DIST_R.fig', ...
               'PLOT/SUB12/SW_DIST_R.fig', ...
               'PLOT/SUB13/SW_DIST_R.fig', ...
               'PLOT/SUB14/SW_DIST_R.fig'};

% Aprire la prima figura per riutilizzare la stessa come base
baseFig = openfig(figurePaths{1}, 'reuse');
baseAxes = findobj(baseFig, 'Type', 'axes');

% Definire i colori per distinguere i segnali provenienti da diverse figure
colors = lines(numel(figurePaths));

% Etichetta la prima figura correttamente
for k = 1:length(baseAxes)
    children1 = get(baseAxes(k), 'Children');
    for j = 1:length(children1)
        if strcmp(get(children1(j), 'LineStyle'), '-') % Solo linee solide
            set(children1(j), 'Marker', 'none'); % Rimuove i marcatori
            set(children1(j), 'LineWidth', 1.5); % Aumenta lo spessore delle linee di SUB1
            set(children1(j), 'Color', colors(1,:)); % Imposta un colore distintivo per SUB1
            set(children1(j), 'DisplayName', 'sub1'); % Imposta il nome per la legenda di SUB1
        end
    end
end

% Loop su tutte le altre figure
for i = 2:numel(figurePaths)
    % Apri la figura successiva
    nextFig = openfig(figurePaths{i}, 'reuse');
    nextAxes = findobj(nextFig, 'Type', 'axes');
    
    % Verifica che il numero di subplot sia uguale
    if length(baseAxes) ~= length(nextAxes)
        error('Le figure devono avere lo stesso numero di subplot.');
    end
    
    % Loop su ogni subplot per combinare i dati
    for k = 1:length(baseAxes)
        % Identifica l'asse corrente in entrambe le figure
        ax1 = baseAxes(k);
        ax2 = nextAxes(k);
        
        % Copia tutti i figli (elementi del plot) da ax2 a ax1
        children2 = get(ax2, 'Children');
        for j = 1:length(children2)
            h = copyobj(children2(j), ax1);
            if strcmp(get(h, 'LineStyle'), '-') % Solo linee solide
                set(h, 'Color', colors(i,:)); % Imposta un colore distintivo per ogni figura
                set(h, 'Marker', 'none'); % Rimuove i marcatori
                set(h, 'LineWidth', 1); % Aumenta lo spessore della linea
                set(h, 'DisplayName', ['sub', num2str(i)]); % Imposta il nome distintivo per la legenda
            end
        end
% Disabilita la griglia per ogni subplot
    grid(baseAxes(k), 'off');
    end
    
    % Chiudere la figura originale dopo aver copiato i contenuti
    close(nextFig);
end

% Aggiungere etichette, titolo e legenda per ogni subplot
for k = 1:length(baseAxes)
    legend(baseAxes(k), 'show');
    xlabel(baseAxes(k), ' tSCS Frequency (Hz)');
    ylabel(baseAxes(k), 'Normalized muscle response(%)');
end

% Salva la nuova figura combinata nella cartella "DATI"
savefig(baseFig, fullfile('SWING', 'SW_DIST_R.fig'));

% Chiudere la figura base
close(baseFig);

