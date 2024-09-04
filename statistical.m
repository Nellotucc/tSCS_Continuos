
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
frequencies = {'0Hz', '5Hz', '15Hz', '30Hz', '50Hz','65', '75Hz', '100Hz'};

%Global group
group = [repmat(frequencies(1), 1, length(group8)), ...
         repmat(frequencies(2), 1, length(group5)), ...
         repmat(frequencies(3), 1, length(group2)), ...
         repmat(frequencies(4), 1, length(group3)), ...
         repmat(frequencies(5), 1, length(group4)), ...
         repmat(frequencies(6), 1, length(group6)), ...
         repmat(frequencies(7), 1, length(group7)), ...
         repmat(frequencies(8), 1, length(group1))];
      

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