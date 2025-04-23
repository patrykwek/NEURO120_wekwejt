% Per-Unit Functional Decoding (Cue vs Lever vs Lick)

% Load behavioral and neural data
cd /Volumes/WD_BLACK
load('ratBEHstruct_unit_firstsection_trial3.mat');
load('celltype_all_sessions.mat');  % Contains: cell_type_all

% Parameters
win = [-500, 2000];       % ms
binSize = 50;             % ms per bin
timeBins = win(1):binSize:win(2);
numBins = length(timeBins) - 1;

% Threshold for minimum spikes per condition
minSpikes = 5;

% Initialize storage
unit_accuracies = [];
unit_labels = {};
confusion_matrices = {};
accuracies_by_event = struct('cue', [], 'lever', [], 'lick', []);
unit_types = [];

% Iterate over all sessions
for sess = 1:length(ratBEHstruct_unit)
    session = ratBEHstruct_unit(sess);
    if ~isfield(session, 'spikes') || isempty(session.spikes), continue; end

    spikes = session.spikes;
    cueTimes = vertcat(session.cuedTimes{:});
    pokeTimes = vertcat(session.pokeTimes{:});
    lickTimes = session.LickOn;

    % Iterate over all channels and units
    for ch = 1:size(spikes, 2)
        units = spikes{1, ch};
        if isempty(units), continue; end

        for u = 1:length(units)
            st = units{u};
            if isempty(st), continue; end

            % Only include classified cells
            if sess > length(cell_type_all) || ch > length(cell_type_all{sess}) || ...
               u > length(cell_type_all{sess}{ch}) || isempty(cell_type_all{sess}{ch}{u}) || cell_type_all{sess}{ch}{u} == 0
                continue;
            end

            this_type = cell_type_all{sess}{ch}{u};

            % Align to each event
            extract_binned = @(times) cell2mat(arrayfun(@(t) histcounts(st - t, timeBins), times, 'UniformOutput', false));
            cue_mat   = extract_binned(cueTimes);
            lever_mat = extract_binned(pokeTimes);
            lick_mat  = extract_binned(lickTimes);

            % Apply spike count threshold
            if sum(cue_mat, 'all') < minSpikes || sum(lever_mat, 'all') < minSpikes || sum(lick_mat, 'all') < minSpikes
                continue;
            end

            % Equalize trial counts
            nCue = size(cue_mat, 1);
            nLever = size(lever_mat, 1);
            nLick = size(lick_mat, 1);
            nMin = min([nCue, nLever, nLick]);
            if nMin < 10, continue; end

            cue_mat = cue_mat(1:nMin, :);
            lever_mat = lever_mat(1:nMin, :);
            lick_mat = lick_mat(1:nMin, :);

            % Concatenate and create labels
            X = [cue_mat; lever_mat; lick_mat];
            Y = [zeros(nMin,1); ones(nMin,1); 2*ones(nMin,1)];

            % Train/test split
            cv = cvpartition(length(Y), 'HoldOut', 0.2);
            Xtrain = X(training(cv), :);
            Ytrain = Y(training(cv));
            Xtest = X(test(cv), :);
            Ytest = Y(test(cv));

            % Print label distribution
            fprintf('Unit: Sess%d_Ch%d_U%d\n', sess, ch, u);
            disp('Training label distribution:');
            tabulate(Ytrain)
            disp('Testing label distribution:');
            tabulate(Ytest)

            model = fitcecoc(Xtrain, Ytrain);
            Ypred = predict(model, Xtest);
            acc = mean(Ypred == Ytest);

            % Store
            unit_accuracies(end+1) = acc;
            unit_labels{end+1} = sprintf('Sess%d_Ch%d_U%d', sess, ch, u);
            confusion_matrices{end+1} = confusionmat(Ytest, Ypred);
            unit_types(end+1) = this_type;

            % Store % correct by class
            for class_id = 0:2
                idx = Ytest == class_id;
                if sum(idx) > 0
                    correct = sum(Ypred(idx) == class_id) / sum(idx);
                else
                    correct = NaN;
                end
                switch class_id
                    case 0, accuracies_by_event.cue(end+1) = correct;
                    case 1, accuracies_by_event.lever(end+1) = correct;
                    case 2, accuracies_by_event.lick(end+1) = correct;
                end
            end
        end
    end
end

% Plot distribution of decoding accuracies
figure;
histogram(unit_accuracies, 30);
xlabel('Decoding Accuracy');
ylabel('# Units');
title('Per-Unit Functional Decoding Accuracy (Cue vs Lever vs Lick)');

% Display summary
fprintf('Mean decoding accuracy: %.2f%%\n', mean(unit_accuracies)*100);
fprintf('Units decoded: %d\n', length(unit_accuracies));

% Show first 10 confusion matrices
for i = 1:min(10, length(confusion_matrices))
    figure;
    confusionchart(confusion_matrices{i}, {'Cue','Lever','Lick'});
    title(sprintf('Confusion Matrix - Unit %d (%s)', i, unit_labels{i}));
end

% Plot stacked histograms by event with colors by type
type_colors = [0.2 0.4 1; 0.2 1 0.4; 1 0.4 0.2];  % MSN, FSI, TAN
fields = {'cue', 'lever', 'lick'};
labels = {'Cue', 'Lever Press', 'Lick'};

figure;
bins = 0:5:100;
for i = 1:3
    subplot(1,3,i); hold on;
    accs = accuracies_by_event.(fields{i}) * 100;
    t = unit_types;
    valid = ~isnan(accs);
    accs = accs(valid);
    t = t(valid);

    counts = zeros(3, length(bins)-1);
    for j = 1:3
        acc_by_type = accs(t == j);
        counts(j, :) = histcounts(acc_by_type, bins);
    end

    b = bar(bins(1:end-1)+2.5, counts', 'stacked');
    for j = 1:3
        b(j).FaceColor = type_colors(j,:);
        b(j).EdgeColor = 'k';
    end

    xlabel('% Correct');
    ylabel('# Units');
    title(sprintf('%s Decoding', labels{i}));
    legend({'MSN', 'FSI', 'TAN'});
end
