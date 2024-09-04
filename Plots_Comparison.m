% % clear all
% % close all
% % clc
% % 
% % 
% % % Percorsi delle figure
% % path1 = 'SWING/SW_DIST_L_GROUP.fig';
% % path2 = 'SWING/SW_DIST_R_GROUP.fig';
% % 
% % 
% % % Aprire la prima figura
% % fig1 = openfig(path1, 'reuse'); % 'reuse' per riutilizzare la stessa figura senza creare nuove finestre
% % fig2 = openfig(path2, 'reuse');
% % 
% % 
% % % Trova gli oggetti grafici in entrambe le figure
% % axes1 = findobj(fig1, 'Type', 'axes');
% % axes2 = findobj(fig2, 'Type', 'axes');
% % 
% % 
% % % Assicurati che il numero di assi in entrambe le figure sia lo stesso
% % if length(axes1) ~= length(axes2)
% %     error('Le due figure devono avere lo stesso numero di subplot');
% % end
% % 
% % % Definisci i colori per i segnali
% % colors = lines(2); % Usare la funzione lines per colori distinti
% % 
% % % Copiare gli oggetti grafici da ogni subplot di fig2 a fig1
% % for k = 1:length(axes2)
% %     % Identifica l'asse corrente in entrambe le figure
% %     ax1 = axes1(k);
% %     ax2 = axes2(k);
% %     
% %     % Aggiorna il colore dei segnali della prima figura
% %     children1 = get(ax1, 'Children');
% %     for j = 1:length(children1)
% %         if strcmp(get(children1(j), 'LineStyle'), '-') % Solo linee solide
% %             set(children1(j), 'Color', colors(1,:)); % Imposta il colore per i segnali della prima figura
% %             set(children1(j), 'DisplayName', 'L'); % Imposta il nome per la legenda
% %         end
% %     end
% %     
% %     % Copia tutti i figli (elementi del plot) da ax2 a ax1
% %     children2 = get(ax2, 'Children');
% %     for j = 1:length(children2)
% %         h = copyobj(children2(j), ax1);
% %         if strcmp(get(h, 'LineStyle'), '-') % Solo linee solide
% %             set(h, 'Color', colors(2,:)); % Imposta il colore per i segnali della seconda figura
% %             set(h, 'DisplayName', 'R'); % Imposta il nome per la legenda
% %         end
% %     end
% % end
% % 
% % % Aggiungere etichette, titolo e legenda se necessario
% % % Nota: Questo può essere fatto individualmente per ogni asse se necessario
% % for k = 1:length(axes1)
% %     legend(axes1(k), 'show');
% %     xlabel(axes1(k), 'Frequencies (Hz)');
% %     ylabel(axes1(k), 'Normalized muscle activity (%)');
% %     title(axes1(k), ['Combined Plot ', num2str(k)]);
% % end
% % 
% % % Salva la nuova figura combinata
% % savefig(fig1, 'combined_figure.fig');
% % 
% % % Chiudere le figure originali se non servono più
% % close(fig2);
% % close(fig1);



clear all
close all
clc

% Percorsi delle figure
path1 = 'STANCE/ST_DIST_L_GROUP.fig';
path2 = 'STANCE/ST_DIST_R_GROUP.fig';

% Apri le figure
fig1 = openfig(path1, 'reuse');
fig2 = openfig(path2, 'reuse');

% Trova gli oggetti grafici in entrambe le figure
axes1 = findobj(fig1, 'Type', 'axes');
axes2 = findobj(fig2, 'Type', 'axes');

% Assicurati che il numero di assi in entrambe le figure sia lo stesso
if length(axes1) ~= length(axes2)
    error('Le due figure devono avere lo stesso numero di subplot');
end

% Definisci i colori per i segnali
colors = lines(2); % Usare la funzione lines per colori distinti

% Definisci l'ordine desiderato per combinare i subplot
% Esempio: voglio combinare il 1° subplot di fig1 con il 3° di fig2, e così via.
subplot_pairs = [1, 7; 2, 1; 3, 2; 4, 3; 5, 4; 6, 5; 7, 6];  % [indice_fig1, indice_fig2]

% Copia e combina i subplot secondo l'ordine specificato
for k = 1:size(subplot_pairs, 1)
    % Indici dei subplot da combinare
    idx1 = subplot_pairs(k, 1);
    idx2 = subplot_pairs(k, 2);
    
    % Identifica l'asse corrente in entrambe le figure
    ax1 = axes1(idx1);
    ax2 = axes2(idx2);
    
    % Aggiorna il colore dei segnali della prima figura
    children1 = get(ax1, 'Children');
    for j = 1:length(children1)
        if strcmp(get(children1(j), 'LineStyle'), '-') % Solo linee solide
            set(children1(j), 'Color', colors(1,:)); % Imposta il colore per i segnali della prima figura
            set(children1(j), 'DisplayName', 'L'); % Imposta il nome per la legenda
        end
    end
    
    % Copia tutti i figli (elementi del plot) da ax2 a ax1
    children2 = get(ax2, 'Children');
    for j = 1:length(children2)
        h = copyobj(children2(j), ax1);
        if strcmp(get(h, 'LineStyle'), '-') % Solo linee solide
            set(h, 'Color', colors(2,:)); % Imposta il colore per i segnali della seconda figura
            set(h, 'DisplayName', 'R'); % Imposta il nome per la legenda
        end
    end
end

% Aggiungere etichette, titolo e legenda se necessario
for k = 1:size(subplot_pairs, 1)
    idx1 = subplot_pairs(k, 1);
    legend(axes1(idx1), 'show');
    xlabel(axes1(idx1), 'Frequencies (Hz)');
    ylabel(axes1(idx1), 'Normalized muscle activity (%)');
    title(axes1(idx1), ['Combined Plot ', num2str(k)]);
end

% Salva la nuova figura combinata
savefig(fig1, 'combined_figure.fig');

% Chiudi le figure originali se non servono più
close(fig2);
close(fig1);

