%% Clean-up
clear; 
close all;
% warning('off', 'MATLAB:print:FigureTooLargeForPage')
% Specify parameters
Options.directory = 'D:\JG\Project Gi-DREADD\Auswertung\P11_Emx1_NaCl\';
Options.lengthFOVid = 12; % number of characters specifying animal & FOV
Options.timestamp = datestr(now, 'yymmdd_HHMMSS');
Options.age = 'P11'; % specify age group
Options.mouseline = 'Emx1'; % specify mouse line
Options.subtance = 'NaCl'; % specify substance that was injected
Options.based_on = mfilename('fullpath');
cond2 = {'control', Options.subtance};
timewindow.control = [0 30]; % time window to calculate average values for statistics
timewindow.C21 = [30 60]; % time window to calculate average values for statistics
SPW_timewindow.seconds = 300; % length of time window for SPW occurence frequency progression 
SPW_timewindow.control = 30; % total length for control
SPW_timewindow.C21 = 120; % total length for C21 
newbandpower = [8 40];
% Get file list
Options.FileList = dir([Options.directory, '*JG', '*.mat']);
N_animals = size(Options.FileList, 1);
Documentation = cell(6, N_animals); % number of parameters
Analysis.SPWprogression = NaN(N_animals,(SPW_timewindow.control/SPW_timewindow.seconds*60+SPW_timewindow.C21/SPW_timewindow.seconds*60)); %prepare matrix with rows: number of animals; columns: number of time windows defined by parameters in SPW_timewindow
for x = 1:N_animals
    Documentation{1, x} = Options.FileList(x).name(1:Options.lengthFOVid);
end
% check for duplicated FOV IDs
aux = strjoin(Documentation(1, :)); 
for x = 1:N_animals
    if length(strfind(aux, Documentation{1, x})) > 1
        error(['Error: The animal with running #', num2str(x), ' appears more than once.'])
    end
end
clear aux
% add animal-related file names
for x = 1:N_animals
    Documentation{2, x}.Files = Options.FileList(x);
    Documentation{3, x} = load([Options.directory, Documentation{2, x}.Files.name], ...
        'options', 'metadata', 'bands', 'SPW','bp'); %load these specified variables for documentation
    
end

disp([num2str(N_animals), ' animals in that group were found ...'])
% animal-wise computations
tic
for x = 1:N_animals
    animalID = Documentation{1, x};
    disp(['animal # ', num2str(x), ' is being processed ...'])
    load([Options.directory, Documentation{2, x}.Files.name]); % loading the .mat file from that animal

    % ========== Counts and frequencies ==========
    Options.SPW_thr(x) = SPW_thr;
    Analysis.AnimalID{x,1} = animalID;
    for n = 3 % for LFP_Ref
        field = string(cond{n});
        bands.(field){7,1} = newbandpower; 
        FinalTime.(field)(x,1) = metadata.rectime{2, n}; % in s; after drift removal
        DriftPercent.(field)(x,1) = (metadata.rectime{1, n}-metadata.rectime{2, n})/metadata.rectime{1, n}; % Amount of Drift in % that was removed
        Analysis.pxx.(field){1,1}(x,:) = pxx{n};
        for i = 1:size(bands.(field),1)            
            bp.LFP_Ref(i,1) = bandpower(pxx{n}, w, bands.(field){i,1} , 'psd'); % Bandpower for all specified bands
            Analysis.bp_raw{i,1}(x,n) = bp.LFP_Ref(i,1);
        end
    end
    for n = 1:(length(cond)-1) % for control and C21 without LFP_Ref
        field = string(cond{n});
        FinalTime.(field)(x,1) = metadata.rectime{2, n}; % in s; after drift removal
        DriftPercent.(field)(x,1) = (metadata.rectime{1, n}-metadata.rectime{2, n})/metadata.rectime{1, n}; % Amount of Drift in % that was removed
        Analysis.SPW(x,n) = length(SPW.(field){4, 1})/FinalTime.(field)(x,1)*60; % SPW occurence in 1/min
        Analysis.pxx.(field){1,1}(x,:) = pxx{n}; % PSD of entire recording
        Analysis.w = w';
        Analysis.osci_w = osci.w{1,1}';
        % PSD of subsections of the recording to see if power increases monotonically during recording
        for o = 1:subsec.n{n}
            subsec.bp.(field){o} = bandpower(subsec.pxx.(field){o}, w, newbandpower , 'psd'); % recalculate bandpower for the new frequency range
            Analysis.subsec_bp.(field)(x,o) = subsec.bp.(field){o}; % extract bandpower of subsections of length subsec.length
        end
        % PSD of specific time windows for control and C21
        timewindow.datapoints.(field) = timewindow.(field)*60*metadata.lfp.fs;
        win_begin = timewindow.datapoints.(field)(1)+1;
        win_end = timewindow.datapoints.(field)(2);
        if win_end > length(lfp{1,n})
            win_end = length(lfp{1,n});
        end
        win_length = win_welch*metadata.lfp.fs; % larger windows increase resolution for smaller frequency ranges
        [pxx{n}, w] = pwelch(lfp{1,n}(win_begin:win_end), win_length, noverlap,[], metadata.lfp.fs, 'psd');   
        bands.(field){7,1} = newbandpower;     
        for i = 1:size(bands.(field),1)
            Analysis.bp_raw{i,1}(x,n) = bandpower(pxx{n}, w, bands.(field){i,1} , 'psd'); % Bandpower for all specified bands
            Analysis.bp_diff{i,1}(x,n) = Analysis.bp_raw{i,1}(x,n)-Analysis.bp_raw{i,1}(x,3); % bandpower difference between condition and LFP_Ref
            Analysis.bp_normalized{i,1}(x,n) = Analysis.bp_raw{i,1}(x,n)/Analysis.bp_raw{i,1}(x,3); % bandpower for each condition normalized to LFP_Ref
        end
    end
    Documentation{3, x}.bands = bands;
% SPW occurrence frequency for specific time windows    
    for n = 1:(length(cond)-1) % for control and C21 only
        field = string(cond{n});
        win_begin = timewindow.datapoints.(field)(1)+1;
        win_end = timewindow.datapoints.(field)(2);        
        n_SPW = SPW.(field){4, 1}(win_begin<=SPW.(field){4, 1} & SPW.(field){4, 1}<win_end);
        Analysis.SPW_timewindow(x,n) = length(n_SPW)/diff(timewindow.(field));
    end

% LFP burst occurrence frequency for specific time windows    
    for n = 1:(length(cond)-1) % for control and C21 only
        field = string(cond{n});
        timewindow.seconds.(field) = timewindow.(field)*60;
        win_begin = timewindow.seconds.(field)(1)+1;
        win_end = timewindow.seconds.(field)(2);
        for i = 1:size(bands.(field),1)-1 %only bands 1-6
            n_LFPbursts_timewindow = bands.(field){i, 16}(win_begin<=bands.(field){i, 16} & bands.(field){i, 16}<win_end);
            Analysis.LFPbursts_timewindow{i,1}(x,n) = length(n_LFPbursts_timewindow)/diff(timewindow.(field));
        end
    end

% SPW occurrence frequency progression 
    for n = 1:(length(cond)-1) % for control and C21 only
        field = string(cond{n});
        datapoints = SPW_timewindow.seconds*metadata.lfp.fs; % datapoints of SPW_timewindow, for extracting n_SPW
        n_windows.(field) = min([floor(metadata.rectime{2, n}/SPW_timewindow.seconds), SPW_timewindow.(field)/SPW_timewindow.seconds*60]); % calculate n_windows in dependence of SPW_timewindow.seconds, with maximum n_windows given by parameters set in SPW_timewindow
        q.(field) = (1+(n-1)*SPW_timewindow.control/SPW_timewindow.seconds*60):1:(n_windows.(field)+(n-1)*SPW_timewindow.control/SPW_timewindow.seconds*60); % q gives column to write SPW occurence frequency in Analysis.SPWprogression
        for m = 1:n_windows.(field)
            win_begin = 1+datapoints*(m-1);
            win_end = datapoints*m;
            n_SPW = SPW.(field){4, 1}(win_begin<=SPW.(field){4, 1} & SPW.(field){4, 1}<win_end);            
            Analysis.SPWprogression(x,q.(field)(m)) = length(n_SPW)/SPW_timewindow.seconds*60;
        end
    end
end
%
%Prepare SPW progression and cond_SPW for plotting
Analysis.SPWprogression_mean = nanmean(Analysis.SPWprogression, 1);
Analysis.SPWprogression_sem = nanstd(Analysis.SPWprogression,0,1)/sqrt(size(Analysis.SPWprogression,1));
q = 1;
for n = 1:(length(cond)-1) % for control and C21 only
    field = string(cond{n});
    b = SPW_timewindow.(field)/SPW_timewindow.seconds*60; %given by parameters set in SPW_timewindow
    for p = 1:b
        cond_SPW{1,q} = [cond2{n}, '_', num2str(p)];
        q = q+1;
    end
end


for n = 1:length(cond)
% Calculate average PSD across animals
    field = string(cond{n});
    Analysis.pxx.(field){1,2} = mean(Analysis.pxx.(field){1,1}, 1);
    Analysis.pxx.(field){1,3} = std(Analysis.pxx.(field){1,1},0,1)/sqrt(size(Analysis.pxx.(field){1,1},1));
% Calculate average PSD of bandpower across animals
%     for i = 1:size(bands.(field),1)            
%         Analysis.osci_pxx.(field){i,2} = mean(Analysis.osci_pxx.(field){i,1}, 1);
%         Analysis.osci_pxx.(field){i,3} = std(Analysis.osci_pxx.(field){i,1},0,1)/sqrt(size(Analysis.osci_pxx.(field){i,1},1));
%     end
end
% Calculate normalized PSD for each condition to LFP_Ref
for n = 1:length(cond)-1
    field = string(cond{n});
    Analysis.pxx.(field){1,4} = Analysis.pxx.(field){1,1}./Analysis.pxx.LFP_Ref{1,1}; 
    Analysis.pxx.(field){1,5} = mean(Analysis.pxx.(field){1,4}, 1);
    Analysis.pxx.(field){1,6} = std(Analysis.pxx.(field){1,4},0,1)/sqrt(size(Analysis.pxx.(field){1,4},1));
end
% prepare bandpower of subsections
Analysis.subsec_bp.all = [];
q = 1;
for n = 1:(length(cond)-1) % for control and C21 only
    field = string(cond{n});
    Analysis.subsec_bp.(field)(Analysis.subsec_bp.(field)==0) = NaN;
    if n == 1
        b = 6; % fix time for control to 6 subsection (30 min)
    elseif n == 2
        b = 24; % fix time for control to 24 subsection (120 min)
    end
%     b = min(sum(Analysis.subsec_bp.(field)&1, 2));
    Analysis.subsec_bp.all = [Analysis.subsec_bp.all Analysis.subsec_bp.(field)(1:size(Analysis.subsec_bp.(field), 1), 1:b)];
    for p = 1:b
        condsubsec_bp{1,q} = [cond2{n}, '_', num2str(p)];
        q = q+1;
    end
end
Analysis.subsec_bp.all_mean = nanmean(Analysis.subsec_bp.all, 1);
Analysis.subsec_bp.all_sem = nanstd(Analysis.subsec_bp.all,0,1)/sqrt(size(Analysis.subsec_bp.all,1));
Analysis.subsec_bp.LFP_Ref_mean = repelem(mean(Analysis.bp_raw{7, 1}(:,3)), length(Analysis.subsec_bp.all_mean));
Analysis.subsec_bp.LFP_Ref_sem = repelem(std(Analysis.bp_raw{7, 1}(:,3))/sqrt(size(Analysis.bp_raw{7, 1}(:,3), 1)), length(Analysis.subsec_bp.all_mean));

for x = 1:N_animals
    Analysis.subsec_bp.all_normalized(x,:) = Analysis.subsec_bp.all(x,:)/Analysis.bp_raw{7, 1}(x,3);
end
Analysis.subsec_bp.all_normalized_mean = nanmean(Analysis.subsec_bp.all_normalized, 1);
Analysis.subsec_bp.all_normalized_sem = nanstd(Analysis.subsec_bp.all_normalized,0,1)/sqrt(size(Analysis.subsec_bp.all,1));
toc
clear ans bands bp comments driftidx cond_legend drift_ext driftidx driftperiods
clear driftperiods_all eventperiods field goodidx i k x lfp metadata MinEventTime
clear mt_f mt_s mt_s_log mt_t n noverlap nSec_base nSec_win options osci freqs
clear overlap params pxx s Spec_fs SPW SPW_thr std_factor subsec w win win_osci
clear win_osci_f win_overlap win_welch z_mt_s filename file_path
clear datapoints n_windows win_begin win_end n_SPW q o p b 
clear win_length m SPW_thr_STD SPW_thr_STD_factor
% Plots unrelated to different frequency bands
% occurence frequency SPWs
figure
customplot_paired(Analysis.SPW_timewindow, cond2);
ylabel('SPW occurence (1/min)');
title('Occurence of SPWs');
export_fig(([Options.directory, Options.age, '_', Options.mouseline, '_',... 
    Options.subtance, '_', Options.timestamp]), '-pdf', '-silent', '-transparent', '-nocrop'); % '-noinvert'

% progression of occurence frequency SPWs
figure
customplot_paired(Analysis.SPWprogression , cond_SPW);
ylabel('SPW occurence (1/min)');
xtickangle(45)
title('Occurence of SPWs');
export_fig(([Options.directory, Options.age, '_', Options.mouseline, '_',... 
    Options.subtance, '_', Options.timestamp]), '-pdf', '-silent', '-transparent', '-nocrop', '-append');

figure
H = shadedErrorBar([1:1:length(cond_SPW)], Analysis.SPWprogression_mean, Analysis.SPWprogression_sem);
ylabel('SPW occurence (1/min)');
title('Occurence of SPWs');
export_fig(([Options.directory, Options.age, '_', Options.mouseline, '_',... 
    Options.subtance, '_', Options.timestamp]), '-pdf', '-silent', '-nocrop', '-append');

% progression of bandpower of subsections
figure
hold on
customplot_paired(Analysis.subsec_bp.all, condsubsec_bp);
shadedErrorBar([1:1:length(condsubsec_bp)],Analysis.subsec_bp.LFP_Ref_mean, Analysis.subsec_bp.LFP_Ref_mean);
hold off
ylabel('raw Bandpower 8-40 Hz');
set(gca,'YScale','linear');
xtickangle(45)
title('Bandpower for subsections');
export_fig(([Options.directory, Options.age, '_', Options.mouseline, '_',... 
    Options.subtance, '_', Options.timestamp]), '-pdf', '-silent', '-nocrop', '-append');

figure
H = shadedErrorBar([1:1:length(condsubsec_bp)], Analysis.subsec_bp.all_mean, Analysis.subsec_bp.all_sem);
ylabel('raw Bandpower 8-40 Hz');
title('Bandpower for subsections');
export_fig(([Options.directory, Options.age, '_', Options.mouseline, '_',... 
    Options.subtance, '_', Options.timestamp]), '-pdf', '-silent', '-nocrop', '-append');

figure
hold on
H = shadedErrorBar([1:1:length(condsubsec_bp)], Analysis.subsec_bp.all_normalized_mean, Analysis.subsec_bp.all_normalized_sem);
ylabel('normalized Bandpower 8-40 Hz');
title('normalized Bandpower for subsections');
yline(1, '--r')
hold off
export_fig(([Options.directory, Options.age, '_', Options.mouseline, '_',... 
    Options.subtance, '_', Options.timestamp]), '-pdf', '-silent', '-nocrop', '-append');

% Plots for each Frequency band under investigation
% occurance frequency of Events detected in a specific frequency band
figure('Position',[300 100 1000 1000])
tl = tiledlayout('flow','TileSpacing','compact'); 
for i = 1:length(Analysis.LFPbursts_timewindow)            
    nexttile
    hold on
    customplot_paired(Analysis.LFPbursts_timewindow{i,1}, cond2);
    bandlabel = ['[', num2str(Documentation{3, 1}.bands.control{i, 1}(1)),'-',num2str(Documentation{3, 1}.bands.control{i, 1}(2)), ' Hz]'];
    ylabel('Occurence frequency (1/min)', 'Interpreter', 'none');
    title(['Frequency band ', bandlabel]);
    hold off
end
title(tl, [Options.age, '_', Options.mouseline, '_', Options.subtance, ' transient power increase'], 'Interpreter', 'none')
export_fig(([Options.directory, Options.age, '_', Options.mouseline, '_',... 
    Options.subtance, '_', Options.timestamp]), '-pdf', '-silent', '-transparent', '-nocrop', '-append');
%%
% raw Bandpower during control and after injection
figure('Position',[300 100 1000 1000])
tl = tiledlayout('flow','TileSpacing','compact'); 
for i = 1:size(Analysis.bp_raw,1)            
    nexttile
    hold on
    customplot_paired(Analysis.bp_raw{i,1}, [cond2 'LFP_Ref']);
    bandlabel = ['[', num2str(Documentation{3, 1}.bands.control{i, 1}(1)),'-',num2str(Documentation{3, 1}.bands.control{i, 1}(2)), ' Hz]'];
    ylabel('Bandpower', 'Interpreter', 'none');
%     set(gca,'YScale','log');
    title(['Frequency band ', bandlabel]);
    hold off
end
title(tl, [Options.age, '_', Options.mouseline, '_', Options.subtance, ' absolute Bandpower'], 'Interpreter', 'none')
export_fig(([Options.directory, Options.age, '_', Options.mouseline, '_',... 
    Options.subtance, '_', Options.timestamp]), '-pdf', '-silent', '-transparent', '-nocrop', '-append');

% Bandpower during recording minus bandpower during LFP-Ref
figure('Position',[300 100 1000 1000])
tl = tiledlayout('flow','TileSpacing','compact'); 
for i = 1:size(Analysis.bp_raw,1)            
    nexttile
    hold on
    customplot_paired(Analysis.bp_diff{i,1}, cond2);
    bandlabel = ['[', num2str(Documentation{3, 1}.bands.control{i, 1}(1)),'-',num2str(Documentation{3, 1}.bands.control{i, 1}(2)), ' Hz]'];
    ylabel('bp(control)-bp(LFP_Ref)', 'Interpreter', 'none');
    title(['Frequency band ', bandlabel]);
    hold off
end
title(tl, [Options.age, '_', Options.mouseline, '_', Options.subtance, ' Bandpower difference'], 'Interpreter', 'none')
export_fig(([Options.directory, Options.age, '_', Options.mouseline, '_',... 
    Options.subtance, '_', Options.timestamp]), '-pdf', '-silent', '-transparent', '-nocrop', '-append');

% Bandpower during recording normalized bandpower during LFP-Ref
figure('Position',[300 100 1000 1000])
tl = tiledlayout('flow','TileSpacing','compact'); 
for i = 1:size(Analysis.bp_raw,1)            
    nexttile
    hold on
    customplot_paired(Analysis.bp_normalized{i,1}, cond2);
    bandlabel = ['[', num2str(Documentation{3, 1}.bands.control{i, 1}(1)),'-',num2str(Documentation{3, 1}.bands.control{i, 1}(2)), ' Hz]'];
    ylabel('normalized bp', 'Interpreter', 'none');
    title(['Frequency band ', bandlabel]);
    hold off
end
title(tl, [Options.age, '_', Options.mouseline, '_', Options.subtance, ' normalized Bandpower'], 'Interpreter', 'none')
export_fig(([Options.directory, Options.age, '_', Options.mouseline, '_',... 
    Options.subtance, '_', Options.timestamp]), '-pdf', '-silent', '-transparent', '-nocrop', '-append');

% PSD (mean +- SEM) across animals
figure
hold on
for n = 1:length(fieldnames(Analysis.pxx)) % control, C21, LFP_Ref
    field = string(cond{n});
%     if isempty(osci.(field){i,4})
%         continue
%     else
    H = shadedErrorBar(Analysis.w, Analysis.pxx.(field){1,2}, Analysis.pxx.(field){1,3});
    set(H.mainLine, 'Color', cond_color(n,:))
    set(H.mainLine, 'LineWidth', 2)
    set(H.patch, 'FaceColor', cond_color(n,:))
    set(H.edge, 'Color', cond_color(n,:))
end
legend([cond2 'LFP_Ref'], 'Interpreter', 'none')
xlabel('Frequency [Hz]')
ylabel('Power')
set(gca,'YScale','log'); % set(gca,'YScale','linear');
xlim([3 40])
title('average PSD using pwelch')
hold off
export_fig(([Options.directory, Options.age, '_', Options.mouseline, '_',... 
    Options.subtance, '_', Options.timestamp]), '-pdf', '-silent', '-nocrop', '-append'); %'-noinvert'

% normalized PSD (mean +- SEM) across animals
figure
hold on
for n = 1:length(fieldnames(Analysis.pxx))-1 % control, C21
    field = string(cond{n});
    H = shadedErrorBar(Analysis.w, Analysis.pxx.(field){1,5}, Analysis.pxx.(field){1,6});
    set(H.mainLine, 'Color', cond_color(n,:))
    set(H.mainLine, 'LineWidth', 2)
    set(H.patch, 'FaceColor', cond_color(n,:))
    set(H.edge, 'Color', cond_color(n,:))
end
yline(1, '--r')
legend([cond2], 'Interpreter', 'none')
xlabel('Frequency [Hz]')
ylabel('normalised Power')
% set(gca,'YScale','log'); 
set(gca,'YScale','linear');
xlim([3 40])
title('average normalized PSD using pwelch')
hold off
export_fig(([Options.directory, Options.age, '_', Options.mouseline, '_',... 
    Options.subtance, '_', Options.timestamp]), '-pdf', '-silent', '-nocrop', '-append'); %'-noinvert'
% print(gcf, '-dpdf', 'transparent.pdf');

% saving
clear a b H i n tl field bandlabel animalID
save([Options.directory, Options.age, '_', Options.mouseline, '_',... 
    Options.subtance, '_', Options.timestamp, '.mat'], '-v7.3', '-nocompression');
disp(['File was saved to ', ...
    Options.directory, Options.age, '_', Options.mouseline, '_',... 
    Options.subtance, '_', Options.timestamp, '.mat'])
