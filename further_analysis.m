% Load z-scored data and classification
load('zFR_all_sessions_aligned.mat');  % zFR_cue_all, zFR_lick_all, zFR_lever_all
load('celltype_all_sessions.mat');  % cell_type_all

% Build neuron classification vector aligned with zFR matrices
neuronTypes_all = [];
zCue_all = zFR_cue_all;
zLick_all = zFR_lick_all;
zLever_all = zFR_lever_all;

% Only add neurons that were actually included in zFR matrices
count = 0;
for sess = 1:length(cell_type_all)
    for ch = 1:length(cell_type_all{sess})
        for u = 1:length(cell_type_all{sess}{ch})
            val = cell_type_all{sess}{ch}{u};
            if ~isempty(val) && val > 0
                count = count + 1;
                if count <= size(zCue_all, 1)
                    neuronTypes_all(count) = val;
                end
            end
        end
    end
end

% Equalize neuron counts per condition
minN = min([size(zCue_all,1), size(zLick_all,1), size(zLever_all,1)]);
zCue   = zCue_all(1:minN, :);
zLick  = zLick_all(1:minN, :);
zLever = zLever_all(1:minN, :);

% Concatenate features and create labels
X = [zCue; zLick; zLever];
Y = [zeros(minN,1); ones(minN,1); 2*ones(minN,1)];

% Train/test split
cv = cvpartition(length(Y), 'HoldOut', 0.3);
Xtrain = X(training(cv), :);
Ytrain = Y(training(cv));
Xtest = X(test(cv), :);
Ytest = Y(test(cv));

% Train multiclass SVM
model = fitcecoc(Xtrain, Ytrain);
Ypred = predict(model, Xtest);

% Accuracy and confusion matrix
accuracy = mean(Ypred == Ytest);
fprintf('Decoding accuracy: %.2f%%\n', accuracy * 100);
confMat = confusionmat(Ytest, Ypred);
disp('Confusion matrix:');
disp(confMat);

% Plot confusion matrix
fig1 = figure;
confusionchart(confMat, {'Cue','Lick','Lever'});
title('Population Decoding Confusion Matrix');
saveas(fig1, '/Volumes/WD_BLACK/NEURO120/plots/confusion_matrix.png');

%% --- Sliding Window Decoding ---
windowSize = 8;
stepSize = 1;
binSize = 50;
totalBins = size(zCue, 2);
startBins = 1:stepSize:(totalBins - windowSize + 1);
accuracyOverTime = zeros(length(startBins), 1);

for i = 1:length(startBins)
    idx = startBins(i):(startBins(i)+windowSize-1);
    xCue = zCue(:, idx);
    xLick = zLick(:, idx);
    xLever = zLever(:, idx);
    Xwin = [xCue; xLick; xLever];
    Ywin = [zeros(minN,1); ones(minN,1); 2*ones(minN,1)];

    cv = cvpartition(length(Ywin), 'HoldOut', 0.3);
    Xtrain = Xwin(training(cv), :);
    Ytrain = Ywin(training(cv));
    Xtest = Xwin(test(cv), :);
    Ytest = Ywin(test(cv));

    model = fitcecoc(Xtrain, Ytrain);
    Ypred = predict(model, Xtest);
    accuracyOverTime(i) = mean(Ypred == Ytest);
end

fig2 = figure;
timeMidpoints = (startBins + windowSize/2 - 1) * binSize / 1000;
plot(timeMidpoints, accuracyOverTime * 100, 'k-o', 'LineWidth', 2);
xlabel('Time from Event (s)');
ylabel('Decoding Accuracy (%)');
title('Sliding Window Population Decoding Accuracy');
yline(33.3, 'r--', 'Chance', 'LabelHorizontalAlignment','left');
grid on;
saveas(fig2, '/Volumes/WD_BLACK/NEURO120/plots/sliding_window_decoding.png');

%% --- Per Cell Type Decoding ---
types = {'MSN', 'FSI', 'TAN'};
type_ids = [1, 2, 3];
accuracies = zeros(1, length(type_ids));

zCue = zCue_all;
zLick = zLick_all;
zLever = zLever_all;

for t = 1:length(type_ids)
    thisType = type_ids(t);
    idx = find(neuronTypes_all(1:minN) == thisType);

    if length(idx) < 5
        fprintf('Skipping %s: not enough classified neurons.\n', types{t});
        continue;
    end

    zCueT = zCue(idx, :);
    zLickT = zLick(idx, :);
    zLeverT = zLever(idx, :);

    minT = min([size(zCueT,1), size(zLickT,1), size(zLeverT,1)]);
    zCueT = zCueT(1:minT,:);
    zLickT = zLickT(1:minT,:);
    zLeverT = zLeverT(1:minT,:);

    X = [zCueT; zLickT; zLeverT];
    Y = [zeros(minT,1); ones(minT,1); 2*ones(minT,1)];

    cv = cvpartition(length(Y), 'HoldOut', 0.3);
    Xtrain = X(training(cv), :);
    Ytrain = Y(training(cv));
    Xtest = X(test(cv), :);
    Ytest = Y(test(cv));

    model = fitcecoc(Xtrain, Ytrain);
    Ypred = predict(model, Xtest);

    acc = mean(Ypred == Ytest);
    accuracies(t) = acc;
    fprintf('%s Decoding Accuracy: %.2f%%\n', types{t}, acc * 100);
end

fig3 = figure;
bar(accuracies * 100);
set(gca, 'XTickLabel', types);
ylabel('Decoding Accuracy (%)');
title('Per-Cell-Type Population Decoding Accuracy');
ylim([0 100]);
grid on;
saveas(fig3, '/Volumes/WD_BLACK/NEURO120/plots/per_celltype_accuracy.png');


