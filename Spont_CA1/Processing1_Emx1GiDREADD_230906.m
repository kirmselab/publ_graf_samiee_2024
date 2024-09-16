%% Comments
% intended for analysis of 2PLSM data from the AOD 2PLSM
%% Clean-up
close all;
warning('off', 'MATLAB:print:FigureTooLargeForPage')
%% Specify parameters
% Options.directory = '\\folder\';
Options.directory_lengthFOVid = 16; % number of characters specifying date/animal/FOV in Options.directory
Options.lengthFOVid = 16; % number of characters specifying ID in *_aligned.mat files
Options.timestamp = datestr(now, 'yymmdd_HHMMSS');
% parameters used in detectGDPs.m:
Options.interval = 500; % time window [ms]
Options.threshold_activecells = 0.2; % threshold fraction of active cells []
%% Get file list
Options.FileListAlignment = dir([Options.directory, '*_aligned.mat']);
Nfovs = size(Options.FileListAlignment, 1);
Documentation = cell(6, Nfovs);
for k = 1:Nfovs
    Documentation{1, k} = Options.FileListAlignment(k).name(1:Options.lengthFOVid);
end
% check for duplicated FOV IDs
aux = strjoin(Documentation(1, :)); 
for k = 1:Nfovs
    if length(strfind(aux, Documentation{1, k})) > 1
        error(['Error: The FOV with running #', num2str(k), ' appears more than once.'])
    end
end
clear aux
% add FOV-related file names
fs = NaN(1, Nfovs); % sampling frequency [Hz]
for k = 1:Nfovs
    FOVid = Documentation{1, k};
    Documentation{2, k}.Files.Alignment = Options.FileListAlignment(k);
    Documentation{2, k}.Files.Events = dir([Options.directory, FOVid, '*_events.mat']);
    if size(Documentation{2, k}.Files.Events, 1) ~= 1
        error(['Error: The FOV with running #', num2str(k), ' has not exactly one *_events.mat file.'])
    end
    Documentation{3, k} = load([Options.directory, Documentation{2, k}.Files.Alignment.name], 'options');
    aux = load([Options.directory, Documentation{2, k}.Files.Alignment.name], 'metadata');
    fs(k) = aux.metadata.recordingfrequency_Hz;
    clear aux
end
disp([num2str(Nfovs), ' field(s) of view (FOV) found ...'])
%% FOV-wise computations
for k = 1:Nfovs
    tic
    FOVid = Documentation{1, k};
    Documentation{6, k} = rng;
    disp(['FOV# ', num2str(k), ' is being processed ...'])
    if exist('Analysis', 'var') == 1; clear Analysis; end
    Events = load([Options.directory, Documentation{2, k}.Files.Events.name]);    
    %% ========== Correct for post hoc discarded ROIs ==========
    fprintf('... Correct for post hoc discarded ROIs ... ')
    if isfield(Events, 'removeROIs')
        if isempty(Events.removeROIs) == 0 % % +++++: added
            aux = cell2mat(Events.removeROIs(:, 1))';
            Documentation{4, k} = aux;
            Events.finalROIs = setdiff(Events.finalROIs, aux);
            Events.non_finalROIs = unique(sort([Events.non_finalROIs, aux]));
        end
    end
    Analysis.NfinalROIs = length(Events.finalROIs);
    centroids = nanmean(Events.coord);
    centroids = [centroids(1:2:end)', centroids(2:2:end)'];
    Analysis.DistFinalROIs = NaN(Events.Nrois, Events.Nrois);
    for m = Events.finalROIs
        for n = Events.finalROIs(Events.finalROIs <= m)
            Analysis.DistFinalROIs(m, n) = sqrt((centroids(m, 1) - centroids(n, 1))^2 + (centroids(m, 2) - centroids(n, 2))^2);
        end
    end
    fprintf('completed\n')
    %% ========== Counts and frequencies ==========
    fprintf('... Counts and frequencies ... ')
    Analysis.FinalTime = length(Events.finalFrames); % [frames]
    Analysis.FinalTime_min = Analysis.FinalTime / fs(k) / 60; % [min]
    Analysis.CaT_Counts = cell(1, Events.Nrois);
    Analysis.CaT_Onsets = cell(1, Events.Nrois);
    Analysis.CaT_Frequencies = cell(1, Events.Nrois);
    for m = Events.finalROIs
        Analysis.CaT_Counts{1, m} = nansum(cell2mat(Events.Results_perROI{12, m}));
        Analysis.CaT_Onsets{1, m} = uint32(find(cell2mat(Events.Results_perROI{12, m}) == 1));
        Analysis.CaT_Frequencies{1, m} = Analysis.CaT_Counts{1, m}/Analysis.FinalTime; % [1/frames]
        Analysis.CaT_Frequencies_permin{1, m} = Analysis.CaT_Frequencies{1, m} * fs(k) * 60; % [1/min]
    end
    fprintf('completed\n')
    %% ========== GDP detection ==========
    fprintf('... GDP detection ... ')
    Ncats = cell2mat(Analysis.CaT_Counts);
    aux1 = NaN(max(Ncats, [], 'omitnan'), numel(Events.finalROIs));
    aux2 = aux1;
    count = 0;
    for m = Events.finalROIs
        count = count + 1;
        if Ncats(count) > 0
            aux1(1:Ncats(count), count) = Analysis.CaT_Onsets{1, m};
            aux2(1:Ncats(count), count) = ones(Ncats(count), 1) .* m;
        end
    end
    data = [reshape(aux1, [], 1), reshape(aux2, [], 1)];
    idx = any(isnan(data), 2);
    data(idx, :) = [];
    data(:, 1) = data(:, 1) .* 1000/fs(k); % convert CaT onsets to ms
    clear aux1 aux2 count idx
    [GDPIndex, GDP_Onsets, GDP_Offsets, data_sorted, options] = detectGDPs(data, numel(Events.finalROIs), Options);
    if numel(GDP_Onsets) == 0
        Analysis.GDP_Count = 0;
        Analysis.GDP_Frequency_permin = 0;
        Analysis.GDP_Onsets_Offsets_ms = [];
        Analysis.GDP_Durations_ms = [];
        Analysis.GDP_Sizes = [];
    else
        Analysis.GDP_Count = numel(GDP_Onsets);
        Analysis.GDP_Frequency_permin = Analysis.GDP_Count / Analysis.FinalTime_min;
        Analysis.GDP_Onsets_Offsets_ms = [GDP_Onsets, GDP_Offsets]; % [ms]
        Analysis.GDP_Durations_ms = Analysis.GDP_Onsets_Offsets_ms(:, 2) - Analysis.GDP_Onsets_Offsets_ms(:, 1); % [ms]
        Analysis.GDP_Sizes = NaN(Analysis.GDP_Count, 1);
        for m = 1:Analysis.GDP_Count
            Analysis.GDP_Sizes(m, 1) = numel(unique(data_sorted(find(GDPIndex == m, 1, 'first'):find(GDPIndex == m, 1, 'last'), 2))) ./ numel(Events.finalROIs);
        end
    end
    fprintf('completed\n')
    %% ========== Plot figures ==========
    fprintf('... Plot figures ... ')
    % ************************************************************
    % Raster plot, drift periods, fraction of active cells
    % ************************************************************
    figure(1)
    scale_width = 0.9;
    scale_heigth = 0.6;
    set(gcf,'Units','normalized','position',[0.05,0.2,scale_width,scale_heigth])
    set(gcf,'renderer','Painters');
    ax1 = subplot(2, 1, 1);
    hold on
    for m = Events.finalROIs
        if ~isempty(Analysis.CaT_Onsets{1, m})
           scatter(Analysis.CaT_Onsets{1, m}, find(Events.finalROIs == m) * ones(1, length(Analysis.CaT_Onsets{1, m})), 2, 'blue')
        end
    end
    area(1:length(Events.driftindex), Events.driftindex * Analysis.NfinalROIs, 'EdgeColor', [1 0 0], 'FaceColor', [1 0 0])
    axis([0, Events.Nframes, 0, Analysis.NfinalROIs])
    title('Raster plot and drift periods')
    xlabel('Frame #')
    ylabel('Cell #')
    hold off
    %
    ax2 = subplot(2, 1, 2);
    hold on
    aux = round(Analysis.GDP_Onsets_Offsets_ms ./ (1000/fs(k)));
    GDPindicator = zeros(1, Events.Nframes);
    if Analysis.GDP_Count > 0
       for m = 1:Analysis.GDP_Count
           GDPindicator(aux(m, 1):aux(m, 2)) = 1;
       end
    end
    plot(1:Events.Nframes, GDPindicator)
    axis([0, Events.Nframes, 0, m])
    title('GDPs')
    xlabel('Frame #')
    ylabel('GDP indicator')
    set(gca,'YTick',[])
    ylim([-0.1, 1.1])
    hold off
    %
    linkaxes([ax1, ax2], 'x' );
    g = pan;
    set(g,'Motion','horizontal','Enable','on');
    h = zoom;
    set(h,'Motion','horizontal','Enable','on');
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    fprintf('completed\n')
    %% ========== Complete documentation ==========
    Documentation{5, k}.SavedFiles.Rasterplot = [Options.directory, 'Processing1_', Options.timestamp, FOVid, '_Rasterplot.fig'];
    Documentation{5, k}.SavedFiles.MAT = [Options.directory, 'Processing1_', Options.timestamp, '_', FOVid, '.mat'];
    %% ========== Save file(s) ==========
    saveas(figure(1), [Options.directory, 'Processing1_', Options.timestamp, '_', FOVid, '_Rasterplot.fig'])
    close(figure(1));
    Analysis = orderfields(Analysis);
    save([Options.directory, 'Processing1_', Options.timestamp, '_', FOVid, '.mat'], 'Analysis', 'Documentation', 'Events', 'FOVid', 'Nfovs', 'Options', 'fs', '-v7.3');
    toc
    disp(['... FOV #', num2str(k), ' completed and saved to ', ...
        Options.directory, 'Processing1_', Options.timestamp, '_', FOVid, '.mat'])
end
warning('on', 'MATLAB:print:FigureTooLargeForPage')
disp('All FOVs finished!')