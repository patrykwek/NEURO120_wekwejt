
% Load data
load('ratBEHstruct_unit_firstsection_trial3.mat');
load('celltype_session2.mat');  % should load variable: cell_type{channel, unit}

% Parameters
win = [-500, 2000];         % time window around event (ms)
binSize = 50;                % ms per bin
smoothSigma = 100;           % Gaussian smoothing width (ms)
timeBins = win(1):binSize:win(2);
binCenters = timeBins(1:end-1) + binSize / 2;
sigmaBins = smoothSigma / binSize;
gaussWin = gausswin(round(sigmaBins * 6));
gaussWin = gaussWin / sum(gaussWin);

% Target alignments
alignments = {'cuedTimes', 'pokeTimes', 'LickOn'};
event_labels = {'Cue', 'Lever Press', 'Lick'};

% Loop over each alignment
for alignIdx = 1:length(alignments)
    fieldname = alignments{alignIdx};
    label = event_labels{alignIdx};

    % Containers
    zFR_types = struct('FSI', [], 'MSN', [], 'TAN', []);
    zFR_all = {};
    zFR_mat = [];

    % ---- Process Session 2 Only ----
    sessionIndex = 2;
    session = ratBEHstruct_unit(sessionIndex);
    spikes = session.spikes;

    % Extract alignment times
    if strcmp(fieldname, 'LickOn')
        alignTimes = session.(fieldname);
    else
        alignTimes = vertcat(session.(fieldname){:});
    end

    for ch = 1:size(spikes, 2)
        neuronGroup = spikes{1, ch};
        if isempty(neuronGroup), continue; end

        for u = 1:length(neuronGroup)
            st = neuronGroup{u};
            if isempty(st), continue; end

            alignedAll = [];
            for i = 1:length(alignTimes)
                aligned = st - alignTimes(i);
                aligned = aligned(aligned >= win(1) & aligned <= win(2));
                alignedAll = [alignedAll; aligned];
            end

            % Compute firing rate and z-score
            fr = histcounts(alignedAll, timeBins) / (length(alignTimes) * (binSize / 1000));
            smoothed = conv(fr, gaussWin, 'same');
            if std(smoothed) == 0, continue; end
            z = (smoothed - mean(smoothed)) / std(smoothed);

            zFR_all{end+1} = z;
            zFR_mat = [zFR_mat; z];

            % Get classification
            this_type = 0;
            if ch <= size(cell_type,1) && u <= size(cell_type,2)
                val = cell_type{ch,u};
                if ~isempty(val)
                    this_type = val;
                end
            end

            % Store by cell type
            switch this_type
                case 1
                    zFR_types.MSN(end+1, :) = z;
                case 2
                    zFR_types.FSI(end+1, :) = z;
                case 3
                    zFR_types.TAN(end+1, :) = z;
            end
        end
    end

    % ---- Plot: Average traces by cell type ----
    figure; hold on;
    if ~isempty(zFR_types.MSN)
        plot(binCenters/1000, mean(zFR_types.MSN, 1), 'b', 'DisplayName', 'MSN', 'LineWidth', 2);
    end
    if ~isempty(zFR_types.FSI)
        plot(binCenters/1000, mean(zFR_types.FSI, 1), 'g', 'DisplayName', 'FSI', 'LineWidth', 2);
    end
    if ~isempty(zFR_types.TAN)
        plot(binCenters/1000, mean(zFR_types.TAN, 1), 'r', 'DisplayName', 'TAN', 'LineWidth', 2);
    end
    if ~isempty(zFR_all)
        allMat = cell2mat(zFR_all');
        plot(binCenters/1000, mean(allMat, 1), 'k--', 'DisplayName', 'All Cells', 'LineWidth', 1.5);
    end

    xline(0, 'k--', 'LineWidth', 1.2);
    xlabel(sprintf('Time from %s (s)', label));
    ylabel('Z-scored Firing Rate');
    title(sprintf('%s-Aligned Avg Z-Scored Activity by Cell Type (Session 2)', label));
    legend;
    grid on;

    % ---- Plot: Heatmap ----
    [~, peakSortIdx] = max(zFR_mat, [], 2);
    [~, unitOrder] = sort(peakSortIdx);

    figure;
    imagesc(binCenters / 1000, 1:size(zFR_mat,1), zFR_mat(unitOrder, :));
    xlabel(sprintf('Time from %s (s)', label));
    ylabel('Unit # (Session 2)');
    title(sprintf('%s-Aligned Z-scored Activity Heatmap (Session 2)', label));
    colormap(flipud(gray));
    colorbar;
    caxis([-1 1]);
    xline(0, 'r', 'LineWidth', 1.2);
end
