%code based on lab-wide code written by Kiah

fs = 30000;  % sampling rate in Hz

cell_type_all = {};  % [session][channel][unit]
peak_to_valley_all = {};
hwhm_all = {};
isi_all = {};
fr_all = {};

for sess_idx = 1:length(avgwaveforms)
    waveforms_sess = avgwaveforms{sess_idx};
    avgISI_sess = avgISI{sess_idx};
    avgFR_sess = avgFR{sess_idx};

    cell_type_all{sess_idx} = cell(size(waveforms_sess));
    peak_to_valley_all{sess_idx} = cell(size(waveforms_sess));
    hwhm_all{sess_idx} = cell(size(waveforms_sess));
    isi_all{sess_idx} = cell(size(waveforms_sess));
    fr_all{sess_idx} = cell(size(waveforms_sess));

    for ch = 1:length(waveforms_sess)
        for u = 1:length(waveforms_sess{ch})
            if ~isempty(waveforms_sess{ch}{u}) && ...
               ch <= length(avgISI_sess) && u <= length(avgISI_sess{ch}) && ...
               ~isempty(avgISI_sess{ch}{u}) && ~isempty(avgFR_sess{ch}{u})

                waveform = waveforms_sess{ch}{u};
                isi_vector = avgISI_sess{ch}{u};  % ISI vector
                isi_val = mean(isi_vector) / fs * 1000;  % ms
                isi_all{sess_idx}{ch}{u} = isi_val;

                fr_val = avgFR_sess{ch}{u};  % Hz
                fr_all{sess_idx}{ch}{u} = fr_val;

                % Process waveform
                waveform_mat = nan(4,64);
                for i = 1:4
                    waveform_mat(i,:) = waveform((i-1)*64 + (1:64));
                end
                [amp,ind] = min(min(waveform_mat,[],2));
                biggest_waveform = waveform_mat(ind,:);
                biggest_waveform = interp1(1:64, biggest_waveform, linspace(1,64,1000));
                micros_per_bin = 1e6/fs*64/1000;

                % Peak-to-valley
                [~, peak_bin] = min(biggest_waveform);
                [~, valley_bin] = max(biggest_waveform(peak_bin:end));
                valley_bin = valley_bin + peak_bin;
                ptv = (valley_bin - peak_bin)*micros_per_bin;
                peak_to_valley_all{sess_idx}{ch}{u} = ptv;

                % HWHM
                baseline = mean(biggest_waveform(1:100));
                wv_range = baseline - amp;
                half_min = wv_range./2 + amp;
                [~, cl_ind_min_left] = min(abs(half_min - biggest_waveform(1:peak_bin)));
                [~, cl_index_right] = min(abs(half_min - biggest_waveform(peak_bin+1:end)));
                cl_ind_min_right = cl_index_right + peak_bin;
                hwhm_val = (cl_ind_min_right - cl_ind_min_left)*micros_per_bin;
                hwhm_all{sess_idx}{ch}{u} = hwhm_val;

                % Classification
                if ~isempty(hwhm_val) && ~isempty(ptv)
                    if hwhm_val > 20 && hwhm_val < 200 && ptv <= 700 && isi_val < 250
                        cell_type_all{sess_idx}{ch}{u} = 2; % FSI
                    elseif hwhm_val > 100 && hwhm_val < 450 && ptv > 400 && isi_val < 650
                        cell_type_all{sess_idx}{ch}{u} = 1; % MSN
                    elseif hwhm_val > 100 && ptv > 450 && isi_val > 300
                        cell_type_all{sess_idx}{ch}{u} = 3; % TAN
                    else
                        cell_type_all{sess_idx}{ch}{u} = 0;
                    end
                else
                    cell_type_all{sess_idx}{ch}{u} = 0;
                end
            end
        end
    end
end

% Save all session classifications
save('celltype_all_sessions.mat', 'cell_type_all', 'peak_to_valley_all', 'hwhm_all', 'isi_all', 'fr_all');
fprintf('Saved classifications for all sessions to "celltype_all_sessions.mat"\n');

% Load and aggregate
load('celltype_all_sessions.mat');

count_unclassified = 0; count_msn = 0; count_fsi = 0; count_tan = 0;
ptv_vals = []; hwhm_vals = []; isi_vals = []; fr_vals = [];

fr_msn = []; fr_fsi = []; fr_tan = [];

for sess = 1:length(cell_type_all)
    for ch = 1:length(cell_type_all{sess})
        for u = 1:length(cell_type_all{sess}{ch})
            val = cell_type_all{sess}{ch}{u};

            % Require all 4 features
            if isempty(val) || ...
               isempty(peak_to_valley_all{sess}{ch}{u}) || ...
               isempty(hwhm_all{sess}{ch}{u}) || ...
               isempty(isi_all{sess}{ch}{u}) || ...
               isempty(fr_all{sess}{ch}{u})
                continue;
            end

            ptv = peak_to_valley_all{sess}{ch}{u};
            hwhm = hwhm_all{sess}{ch}{u};
            isi = isi_all{sess}{ch}{u};
            fr  = fr_all{sess}{ch}{u};

            % Count + collect
            switch val
                case 0, count_unclassified = count_unclassified + 1;
                case 1, count_msn = count_msn + 1; fr_msn(end+1) = fr;
                case 2, count_fsi = count_fsi + 1; fr_fsi(end+1) = fr;
                case 3, count_tan = count_tan + 1; fr_tan(end+1) = fr;
            end

            ptv_vals(end+1) = ptv;
            hwhm_vals(end+1) = hwhm;
            isi_vals(end+1) = isi;
            fr_vals(end+1)  = fr;
        end
    end
end

% ---- Plot classification counts ----
figure;
bar([count_unclassified, count_msn, count_fsi, count_tan]);
xticks(1:4); xticklabels({'Unclassified','MSN','FSI','TAN'});
ylabel('Count');
title('Cell Type Classification');

% ---- Plot feature histograms ----
figure;
subplot(2,2,1);
histogram(ptv_vals, 30);
xlabel('Peak-to-Valley Time (\mus)');
ylabel('Count');
title('Distribution of Peak-to-Valley');

subplot(2,2,2);
histogram(hwhm_vals, 30);
xlabel('HWHM (\mus)');
ylabel('Count');
title('Distribution of HWHM');

subplot(2,2,3);
histogram(isi_vals, 30);
xlabel('Average ISI (ms)');
ylabel('Count');
title('Distribution of Avg ISI');

subplot(2,2,4);
histogram(fr_vals, 30);
xlabel('Avg Firing Rate (Hz)');
ylabel('Count');
title('Distribution of Avg Firing Rate');

% ---- Print firing rate summary per class ----
fprintf('\n--- Firing Rate Summary ---\n');
fprintf('MSN: %.2f ± %.2f Hz (n=%d)\n', mean(fr_msn), std(fr_msn), numel(fr_msn));
fprintf('FSI: %.2f ± %.2f Hz (n=%d)\n', mean(fr_fsi), std(fr_fsi), numel(fr_fsi));
fprintf('TAN: %.2f ± %.2f Hz (n=%d)\n', mean(fr_tan), std(fr_tan), numel(fr_tan));


% Directory to save plots
output_dir = '/Volumes/WD_BLACK/NEURO120/plots';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Save all currently open figures
figHandles = findall(groot, 'Type', 'figure');
for i = 1:length(figHandles)
    figure(figHandles(i));
    filename = fullfile(output_dir, sprintf('figure_%02d.png', i));
    saveas(figHandles(i), filename);
end

fprintf('Saved %d figures to: %s\n', length(figHandles), output_dir);
