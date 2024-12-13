clear; clc; close all;

%% Parametri
fs = 100; % Frequenza di campionamento
positions = {'LowerBack','BackPocket','FrontPocket','Hand','ShoulderBag','CoatPocket'};
positionAbbrev = containers.Map({'LowerBack','BackPocket','FrontPocket','Hand','ShoulderBag','CoatPocket'},...
                                {'LB','BP','FP','H','SB','CP'});

basePath = 'cold data';
outputBasePath = 'hot data';

fprintf('Inizio elaborazione...\n');

for subj = 1:25
    subjectID = sprintf('%03d', subj);
    dataFile = fullfile(basePath, subjectID, 'data.mat');

    fprintf('\n------------------------------------\n');
    fprintf('Elaborazione soggetto %s...\n', subjectID);

    if ~exist(dataFile,'file')
        fprintf('Nessun data.mat per soggetto %s, salto.\n', subjectID);
        continue;
    end

    load(dataFile); % Carica 'data'
    fprintf('Caricato data.mat per soggetto %s.\n', subjectID);

    if ~exist('data','var')
        fprintf('Nessuna variabile data nel soggetto %s. Salto.\n', subjectID);
        continue;
    end
    if ~isfield(data,'TimeMeasure1')
        fprintf('Nessun campo TimeMeasure1 nel soggetto %s. Salto.\n', subjectID);
        continue;
    end

    testFields = fieldnames(data.TimeMeasure1);
    testMask = startsWith(testFields, 'Test');
    testFields = testFields(testMask);
    numTests = length(testFields);

    fprintf('Trovati %d test per soggetto %s.\n', numTests, subjectID);

    if numTests == 0
        fprintf('Nessun test presente per soggetto %s, salto.\n', subjectID);
        continue;
    end

    % Assicura cartella di output
    subjOutputPath = fullfile(outputBasePath, subjectID);
    if ~exist(subjOutputPath, 'dir')
        mkdir(subjOutputPath);
        fprintf('Creata directory di output: %s\n', subjOutputPath);
    end

    % Inizializza i contatori per i bout salvati/scartati (solo per l'ultimo test)
    subjectResults = struct();
    for p = 1:length(positions)
        posName = positions{p};
        subjectResults.(posName).saved = 0;
        subjectResults.(posName).discarded = 0;
    end

    %% Prima parte: ricavo soglie dal Test 5
    test5Name = 'Test5';
    if any(strcmp(testFields, test5Name))
        fprintf('Trovato Test5 per soggetto %s. Calcolo soglie...\n', subjectID);

        test5Data = data.TimeMeasure1.(test5Name);
        trialFields = fieldnames(test5Data);
        trialMask = startsWith(trialFields, 'Trial');
        trialFields = trialFields(trialMask);

        var_dyn = struct(); 
        var_stat = struct();
        for p=1:length(positions)
            var_dyn.(positions{p}) = [];
            var_stat.(positions{p}) = [];
        end

        for tIdx = 1:length(trialFields)
            trialName = trialFields{tIdx};
            trialData = test5Data.(trialName);

            if isfield(trialData,'Standards') && ...
               isfield(trialData.Standards,'INDIP') && ...
               isfield(trialData.Standards.INDIP,'ContinuousWalkingPeriod')

                cwp = trialData.Standards.INDIP.ContinuousWalkingPeriod;
                Starts = [cwp.Start];
                Ends = [cwp.End];

                if isempty(Starts) || isempty(Ends) || length(Starts) ~= length(Ends)
                    fprintf("Incongruenza in Start/End nel Test5 %s del soggetto %s.\n", trialName, subjectID);
                    continue;
                end

                if ~isfield(trialData,'SU')
                    fprintf("Nessun SU in Test5 %s del soggetto %s.\n", trialName, subjectID);
                    continue;
                end

                for p = 1:length(positions)
                    posName = positions{p};
                    if ~isfield(trialData.SU, posName)
                        continue;
                    end
                    posData = trialData.SU.(posName);
                    if ~isfield(posData,'Gyr')
                        continue;
                    end

                    GyrData = posData.Gyr;
                    nSamples = size(GyrData,1);

                    start_idxs = max(1, round(Starts * fs) + 1);
                    end_idxs = round(Ends * fs);
                    [start_idxs, order] = sort(start_idxs);
                    end_idxs = end_idxs(order);
                    end_idxs(end_idxs > nSamples) = nSamples;
                    walk_intervals = [start_idxs(:), end_idxs(:)];

                    % varianze dinamiche
                    var_dyn_bout = [];
                    for b = 1:size(walk_intervals,1)
                        s_idx = walk_intervals(b,1);
                        e_idx = walk_intervals(b,2);
                        if s_idx < e_idx && e_idx <= nSamples
                            segGyr = GyrData(s_idx:e_idx,:);
                            gyrMag = sqrt(segGyr(:,1).^2 + segGyr(:,2).^2 + segGyr(:,3).^2);
                            var_gyr = var(gyrMag);
                            var_dyn_bout = [var_dyn_bout; var_gyr];
                        end
                    end

                    % varianze statiche
                    var_stat_bout = [];
                    all_intervals = [];
                    if ~isempty(walk_intervals)
                        if walk_intervals(1,1) > 1
                            all_intervals = [all_intervals; 1, walk_intervals(1,1)-1];
                        end
                        for b = 1:size(walk_intervals,1)-1
                            prev_e = walk_intervals(b,2);
                            next_s = walk_intervals(b+1,1);
                            if (prev_e+1) < (next_s-1)
                                all_intervals = [all_intervals; prev_e+1, next_s-1];
                            end
                        end
                        if walk_intervals(end,2) < nSamples
                            all_intervals = [all_intervals; walk_intervals(end,2)+1, nSamples];
                        end
                    else
                        % Nessun bout dinamico => tutto statico
                        all_intervals = [1, nSamples];
                    end

                    for si = 1:size(all_intervals,1)
                        s_idx = all_intervals(si,1);
                        e_idx = all_intervals(si,2);
                        if s_idx < e_idx && e_idx <= nSamples
                            segGyr = GyrData(s_idx:e_idx,:);
                            gyrMag = sqrt(segGyr(:,1).^2 + segGyr(:,2).^2 + segGyr(:,3).^2);
                            var_gyr = var(gyrMag);
                            var_stat_bout = [var_stat_bout; var_gyr];
                        end
                    end

                    if ~isempty(var_dyn_bout)
                        var_dyn.(posName) = [var_dyn.(posName); mean(var_dyn_bout)];
                    end
                    if ~isempty(var_stat_bout)
                        var_stat.(posName) = [var_stat.(posName); mean(var_stat_bout)];
                    end
                end
            else
                fprintf("Nessun ContinuousWalkingPeriod in Test5 %s del soggetto %s.\n", trialName, subjectID);
            end
        end

        var_dyn_threshold = struct();
        var_stat_threshold = struct();
        for p=1:length(positions)
            posName = positions{p};
            if ~isempty(var_dyn.(posName))
                var_dyn_threshold.(posName) = mean(var_dyn.(posName));
            else
                var_dyn_threshold.(posName) = NaN;
            end
            if ~isempty(var_stat.(posName))
                var_stat_threshold.(posName) = mean(var_stat.(posName));
            else
                var_stat_threshold.(posName) = NaN;
            end
        end

        fprintf('Soglie dal Test5 (media su trial):\n');
        for p=1:length(positions)
            posName = positions{p};
            fprintf('%s: var_dyn_threshold=%f, var_stat_threshold=%f\n', posName, var_dyn_threshold.(posName), var_stat_threshold.(posName));
        end

    else
        fprintf('Nessun Test5 trovato per soggetto %s, impossibile calcolare soglie.\n', subjectID);
        continue;
    end

    %% Applicazione della logica ai test diversi dal Test5
    for idxTest = 1:numTests
        currentTestName = testFields{idxTest};
        if strcmp(currentTestName, 'Test5')
            % Test5 è servito per le soglie, lo ignoriamo ora
            continue;
        end

        currentTestData = data.TimeMeasure1.(currentTestName);
        fprintf('\n--- Gestione test %s per soggetto %s ---\n', currentTestName, subjectID);

        isLastTest = (idxTest == numTests); % Verifica se è l'ultimo test

        if isfield(currentTestData,'Recording') && ...
           isfield(currentTestData.Recording,'Standards') && ...
           isfield(currentTestData.Recording.Standards,'INDIP') && ...
           isfield(currentTestData.Recording.Standards.INDIP,'ContinuousWalkingPeriod')

            cwp = currentTestData.Recording.Standards.INDIP.ContinuousWalkingPeriod;
            Starts = [cwp.Start];
            Ends = [cwp.End];

            if isempty(Starts) || isempty(Ends) || length(Starts) ~= length(Ends)
                fprintf("Incongruenza in Start/End nel test %s del soggetto %s.\n", currentTestName, subjectID);
                continue;
            end

            if ~isfield(currentTestData.Recording,'SU')
                fprintf("Nessun SU in %s del soggetto %s.\n", currentTestName, subjectID);
                continue;
            end

            start_idxs = max(1, round(Starts * fs) + 1);
            end_idxs = round(Ends * fs);
            [start_idxs, order] = sort(start_idxs);
            end_idxs = end_idxs(order);

            walk_intervals = [start_idxs(:), end_idxs(:)];

            for p=1:length(positions)
                posName = positions{p};
                if ~isfield(currentTestData.Recording.SU, posName)
                    fprintf("Nessun segnale %s in %s per soggetto %s.\n", posName, currentTestName, subjectID);
                    continue;
                end
                posData = currentTestData.Recording.SU.(posName);
                if ~isfield(posData,'Acc') || ~isfield(posData,'Gyr')
                    fprintf("Mancano Acc o Gyr in %s in %s per soggetto %s.\n", posName, currentTestName, subjectID);
                    continue;
                end

                AccData = posData.Acc;
                GyrData = posData.Gyr;
                nSamples = size(AccData,1);
                end_idxs(end_idxs > nSamples) = nSamples;

                posAbbr = positionAbbrev(posName);

                for b = 1:size(walk_intervals,1)
                    s_idx = walk_intervals(b,1);
                    e_idx = walk_intervals(b,2);

                    if s_idx < e_idx && e_idx <= nSamples
                        segGyr = GyrData(s_idx:e_idx,:);
                        gyrMag = sqrt(segGyr(:,1).^2 + segGyr(:,2).^2 + segGyr(:,3).^2);
                        var_gyr = var(gyrMag);

                        if isLastTest
                            % Ultimo test: applico logica soglia
                            stat_th = var_stat_threshold.(posName);
                            dyn_th = var_dyn_threshold.(posName);

                            if isnan(stat_th) || isnan(dyn_th)
                                fprintf("Soglia statica o dinamica non disponibile per %s, impossibile discriminare.\n", posName);
                                continue;
                            end

                            % Stampa dei valori prima della decisione di tagliare
                            fprintf('[%s | Bout %d] var_dyn_test5=%f, var_stat_test5=%f, var_bout_corrente=%f\n', ...
                                posName, b, dyn_th, stat_th, var_gyr);

                            if var_gyr > stat_th
                                % Segmento dinamico, salvo
                                segAcc = AccData(s_idx:e_idx,:);
                                segGyr = GyrData(s_idx:e_idx,:);
                                boutID = num2str(b);
                                baseFilename = sprintf('%s-%s-%s-%s', subjectID, currentTestName, boutID, posAbbr);
                                accFilename = fullfile(subjOutputPath, [baseFilename '_Acc.csv']);
                                gyrFilename = fullfile(subjOutputPath, [baseFilename '_Gyr.csv']);

                                writematrix(segAcc, accFilename);
                                writematrix(segGyr, gyrFilename);
                                fprintf('Ultimo test: %s bout %d var_gyr=%f > %f stat_th, SALVO %s, %s\n', ...
                                    posName, b, var_gyr, stat_th, accFilename, gyrFilename);

                                % Incremento contatore salvati
                                subjectResults.(posName).saved = subjectResults.(posName).saved + 1;
                            else
                                fprintf('Ultimo test: %s bout %d var_gyr=%f <= %f stat_th, NON SALVO\n', ...
                                    posName, b, var_gyr, stat_th);
                                % Incremento contatore scartati
                                subjectResults.(posName).discarded = subjectResults.(posName).discarded + 1;
                            end
                        else
                            % Non ultimo test: nessuna soglia, salvo direttamente i dati
                            segAcc = AccData(s_idx:e_idx,:);
                            segGyr = GyrData(s_idx:e_idx,:);
                            boutID = num2str(b);
                            baseFilename = sprintf('%s-%s-%s-%s', subjectID, currentTestName, boutID, posAbbr);
                            accFilename = fullfile(subjOutputPath, [baseFilename '_Acc.csv']);
                            gyrFilename = fullfile(subjOutputPath, [baseFilename '_Gyr.csv']);

                            writematrix(segAcc, accFilename);
                            writematrix(segGyr, gyrFilename);
                            fprintf('%s %s: salvo direttamente senza logica soglia in %s, %s\n', currentTestName, posName, accFilename, gyrFilename);
                        end

                    end
                end
            end

        else
            % Test precedente all'ultimo: ci sono TrialX
            trialFields = fieldnames(currentTestData);
            trialMask = startsWith(trialFields, 'Trial');
            trialFields = trialFields(trialMask);

            if isempty(trialFields)
                fprintf("Nessun Trial per %s del soggetto %s.\n", currentTestName, subjectID);
                continue;
            end

            for tIdx = 1:length(trialFields)
                trialName = trialFields{tIdx};
                trialData = currentTestData.(trialName);

                if isfield(trialData,'Standards') && ...
                   isfield(trialData.Standards,'INDIP') && ...
                   isfield(trialData.Standards.INDIP,'ContinuousWalkingPeriod')

                    cwp = trialData.Standards.INDIP.ContinuousWalkingPeriod;
                    Starts = [cwp.Start];
                    Ends = [cwp.End];

                    if isempty(Starts) || isempty(Ends) || length(Starts) ~= length(Ends)
                        fprintf("Incongruenza in Start/End in %s %s del soggetto %s.\n", currentTestName, trialName, subjectID);
                        continue;
                    end

                    if ~isfield(trialData,'SU')
                        fprintf("Nessun SU in %s %s del soggetto %s.\n", currentTestName, trialName, subjectID);
                        continue;
                    end

                    start_idxs = max(1, round(Starts * fs) + 1);
                    end_idxs = round(Ends * fs);
                    [start_idxs, order] = sort(start_idxs);
                    end_idxs = end_idxs(order);

                    % Se non è ultimo test, salvo direttamente senza logica soglia
                    for p=1:length(positions)
                        posName = positions{p};
                        if ~isfield(trialData.SU, posName)
                            continue;
                        end
                        posData = trialData.SU.(posName);
                        if ~isfield(posData,'Acc') || ~isfield(posData,'Gyr')
                            continue;
                        end
                        AccData = posData.Acc;
                        GyrData = posData.Gyr;
                        nSamples = size(AccData,1);
                        end_idxs(end_idxs > nSamples) = nSamples;

                        walk_intervals = [start_idxs(:), end_idxs(:)];
                        posAbbr = positionAbbrev(posName);

                        for b = 1:size(walk_intervals,1)
                            s_idx = walk_intervals(b,1);
                            e_idx = walk_intervals(b,2);
                            if s_idx < e_idx && e_idx <= nSamples
                                segAcc = AccData(s_idx:e_idx,:);
                                segGyr = GyrData(s_idx:e_idx,:);
                                boutID = num2str(b);
                                baseFilename = sprintf('%s-%s-%s-%s-%s', subjectID, currentTestName, trialName, boutID, posAbbr);

                                accFilename = fullfile(subjOutputPath, [baseFilename '_Acc.csv']);
                                gyrFilename = fullfile(subjOutputPath, [baseFilename '_Gyr.csv']);

                                writematrix(segAcc, accFilename);
                                writematrix(segGyr, gyrFilename);
                                fprintf('%s %s %s: salvo senza soglia in %s, %s\n', currentTestName, trialName, posName, accFilename, gyrFilename);
                            end
                        end
                    end

                else
                    fprintf("Nessun ContinuousWalkingPeriod in %s %s del soggetto %s.\n", currentTestName, trialName, subjectID);
                end
            end

        end
    end

    % Terminato il loop dei test per questo soggetto, stampo la tabella riepilogativa
    posNames = positions(:);
    savedCounts = zeros(length(positions),1);
    discardedCounts = zeros(length(positions),1);
    for p=1:length(positions)
        posName = positions{p};
        savedCounts(p) = subjectResults.(posName).saved;
        discardedCounts(p) = subjectResults.(posName).discarded;
    end

    T = table(posNames, savedCounts, discardedCounts, ...
        'VariableNames', {'Position', 'Saved', 'Discarded'});

    fprintf('\nTabella riassuntiva per soggetto %s (ultimo test):\n', subjectID);
    disp(T);

    % Salvo la tabella in CSV
    summaryFile = fullfile(subjOutputPath, [subjectID '_Summary.csv']);
    writetable(T, summaryFile);
    fprintf('Tabella salvata in %s\n', summaryFile);
end

fprintf('\nElaborazione completata!\n');

