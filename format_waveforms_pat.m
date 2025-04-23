%code based on lab-wide code written by Cheshta Bhatia

cd /Volumes/WD_BLACK/NEURO120/NEURO120_wekwejt
%load('spikeshape_section1.mat');
%load('spikesp17_2.mat');

load('spikeshape_section1.mat');
load('spikesp17_1.mat');

ephys_struct = spikeshape;
Firingrates=cell(1,length(ephys_struct));
spikes2 = cell(1,length(ephys_struct));
spikes3 = cell(1,length(ephys_struct));

for sess = 1:length(ephys_struct)
    avgwaveforms{sess} = cell(1,length(ephys_struct{sess}));
    avgwaveformsv2{sess} = cell(1,length(ephys_struct{sess}));
end
for sess = 1:length(ephys_struct)
    for ch = 1:length(ephys_struct{sess})
        avgwaveforms{sess}{ch} = cell(1,length(ephys_struct{sess}{ch}));
        sizeunits{sess}{ch} = numel(spikeshape{sess}{ch});
    end
end

numUnits = max(cell2mat(sizeunits{1,2}));
numUnits2 = sum(cell2mat(sizeunits{1,2}));
for sess = 1:length(ephys_struct)
    for ch = 1:length(ephys_struct{sess})
        avgwaveforms{sess}{ch} = cell(1,numUnits);
        avgwaveformsv2{sess}{ch} = cell(1,numel(spikeshape{sess}{ch}));
        Firingrates{sess}{ch} = cell(1,numUnits);
    end
end

spikes2 = cell(1,length(avgwaveforms));
for sess = 1:length(ephys_struct)
    for ch = 1:length(ephys_struct{sess})
        spikes2{sess}{ch} = cell(1,numUnits);
    end
end

for sess = 1:length(spikesp17)
    for ch = 1:length(spikesp17{sess})
        for u = 1:length(spikesp17{sess}{ch})
            spikes2{sess}{ch}{u} = spikesp17{sess}{ch}{u};
            spikes3{sess}{ch}{u} = spikesp17{sess}{ch}{u};
        end
    end
end

% Initialize avgISI with the same structure as avgwaveforms
avgISI = cell(size(spikes2));  % {sess}{ch}{unit}

for sess = 1:length(spikes2)
    for ch = 1:length(spikes2{sess})
        for u = 1:length(spikes2{sess}{ch})
            spiketrain = spikes2{sess}{ch}{u};
            if ~isempty(spiketrain) && numel(spiketrain) > 1
                isi = diff(spiketrain);
                avgISI{sess}{ch}{u} = mean(isi);  % Store mean ISI
            else
                avgISI{sess}{ch}{u} = [];
            end
        end
    end
end

avgFR = cell(size(spikes2));  % Match structure of avgISI

for sess = 1:length(spikes2)
    for ch = 1:length(spikes2{sess})
        for u = 1:length(spikes2{sess}{ch})
            spiketrain = spikes2{sess}{ch}{u};
            if ~isempty(spiketrain)
                avgFR{sess}{ch}{u} = length(spiketrain) / 3600;  % Hz
            else
                avgFR{sess}{ch}{u} = [];
            end
        end
    end
end




for sess=1:length(spikes3)
    spikes2v2{sess} = {};
    for i = 1:numel(spikes3{sess})
        spikes2v2{sess} = [spikes2v2{sess}, spikes3{sess}{i}];
    end
end

for sess = 1:length(ephys_struct)
    if ~isempty(ephys_struct{sess})
        for ch = 1:length(ephys_struct{sess})
            if ~isempty(ephys_struct{sess}{ch})
                for unit = 1:length(ephys_struct{sess}{ch})
                    if ~isempty(ephys_struct{sess}{ch}{unit})
                        avgwaveforms{sess}{ch}{unit} = mean(spikeshape{sess}{ch}{unit},2);
                        avgwaveformsv2{sess}{ch}{unit} = mean(spikeshape{sess}{ch}{unit},2);
                    end
                end
            end
        end
    end
end

for sess=1:length(avgwaveformsv2)
    if ~isempty(avgwaveformsv2{sess})
        avgwaveformsv3{sess} = cell(1,numUnits2);
    end
end

for sess=1:length(avgwaveformsv2)
    avgwaveformsv3{sess} = {};
    for i = 1:numel(avgwaveformsv2{sess})
        avgwaveformsv3{sess} = [avgwaveformsv3{sess}, avgwaveformsv2{sess}{i}];
    end
end

numChannels=16;
waveforms = cell(numChannels,numUnits);
avgFR_perChanUnit = cell(numChannels, numUnits);
prop_isi = cell(numChannels,numUnits);

waveformsv2 = cell(sess,numUnits2);
avgFRv2 = cell(sess,numUnits2);
prop_isiv2 = cell(sess,numUnits2);

avgArray = cell(numChannels, numUnits);
for ch = 1:numChannels
    for u = 1:numUnits
        temp = zeros(256,1);count=0;
        for sess = 1:length(avgwaveforms)
            if ~isempty(avgwaveforms{sess})
                if ~isempty(avgwaveforms{sess}{ch}{u})
                    count = count+1;
                    temp = temp + avgwaveforms{sess}{ch}{u};
                    waveforms{ch,u} = temp./count;  
                    waveformsv2{sess,u} = temp./count;  
                end
            end
        end
    end
end

for ch = 1:numChannels
    for u = 1:numUnits
        temp = 0;count=0;
        for sess = 1:length(Firingrates)
            if ~isempty(spikes2{sess})
                if ~isempty(spikes2{sess}{ch}{u})
                    count = count+1;
                    temp = temp + length(spikes2{sess}{ch}{u});
                    avgFR3{ch,u} = temp./3600;
                end
            end
        end
    end
end

for ch = 1:numChannels
    for u = 1:numUnits
        temp = 0;count=0;
        for sess = 1:length(Firingrates)
            if ~isempty(spikes2{sess})
                if ~isempty(spikes2{sess}{ch}{u})
                    count = count+1;
                    spiketrain = spikes2{sess}{ch}{u};
                    isi = diff(spiketrain);
                    prop_isi{ch,u} = sum(isi(isi > 2000))/range(double(spiketrain));
                end
            end
        end
    end
end

avgwaveformsv4 = cell(length(avgwaveformsv3),numUnits2);
for sess = 1:length(avgwaveformsv3)
    if ~isempty(avgwaveformsv3{sess})
        for u = 1:numUnits2
            if ~isempty(avgwaveformsv3{sess}{u})
                avgwaveformsv4{sess,u} = avgwaveformsv3{sess}{u};
            end
        end
    end
end

for u = 1:numUnits2
    temp = 0;count=0;
    for sess = 1:length(Firingrates)
        if ~isempty(spikes2v2{sess})
            if ~isempty(spikes2v2{sess}{u})
                count = count+1;
                temp = temp + length(spikes2v2{sess}{u});
                avgFRv2{sess,u} = temp./3600;
            end
        end
    end
end

for u = 1:numUnits2
    temp = 0;count=0;
    for sess = 1:length(Firingrates)
        if ~isempty(spikes2v2{sess})
            if ~isempty(spikes2v2{sess}{u})
                count = count+1;
                spiketrain = spikes2v2{sess}{u};
                isi = diff(spiketrain);
                prop_isiv2{sess,u} = sum(isi(isi > 2000))/range(double(spiketrain));
            end
        end
    end
end

waveforms_all = {};
avgFRv2_all = {};
prop_isi_all = {};

for sess = 1:size(avgwaveformsv4,1)
    for u = 1:size(avgwaveformsv4,2)
        if ~isempty(avgwaveformsv4{sess,u})
            waveforms_all{end+1} = avgwaveformsv4{sess,u};
            avgFRv2_all{end+1} = avgFRv2{sess,u};
            prop_isi_all{end+1} = prop_isiv2{sess,u};
        end
    end
end

% Directory to save waveform plots
output_dir = '/Volumes/WD_BLACK/NEURO120/waveform_analysis/avg_waveforms';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

plot_count = 1;

for sess = 1:length(avgwaveforms)
    for ch = 1:length(avgwaveforms{sess})
        for unit = 1:length(avgwaveforms{sess}{ch})
            wf = avgwaveforms{sess}{ch}{unit};
            if ~isempty(wf)
                fig = figure('Visible', 'off');  % Do not display the figure
                plot(wf);
                title(sprintf('Session %d, Channel %d, Unit %d', sess, ch, unit));
                xlabel('Sample');
                ylabel('Amplitude');
                axis tight;

                % Build filename and save as PNG
                filename = sprintf('sess%d_ch%d_unit%d.png', sess, ch, unit);
                saveas(fig, fullfile(output_dir, filename));

                close(fig);  % Close the figure to save memory
                plot_count = plot_count + 1;
            end
        end
    end
end

fprintf('Saved %d waveform plots to: %s\n', plot_count - 1, output_dir);


