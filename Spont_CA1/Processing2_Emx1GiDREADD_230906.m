%% Clean-up
clearvars; close all;
%% Specify parameters
Options2.directory = '\\folder\';
Options2.subdirectories = {'folder_groupA\', 'folder_groupB\'}; % cell; each sub-directory represents one group in analysis
Options2.lengthFOVid = 16; % number of characters specifying animal & FOV
Options2.directory_out = [Options2.directory, 'folder_out\']; % root output directory where a new subfolder will be added
%
Options2.edgesF = 0:0.01:20; % bin edges of cumulative distributions of CaT frequency [1/min]
Options2.timestamp = datestr(now, 'yymmdd_HHMMSS');
%% Get file list ('Processing1_*.mat') & make output directory
warning('off', 'MATLAB:print:FigureTooLargeForPage')
Ngroups = length(Options2.subdirectories);
Options2.FileLists = cell(1, Ngroups);
for k = 1:length(Options2.subdirectories)
    Options2.FileLists{k} = dir([Options2.directory, Options2.subdirectories{k}, 'Processing1_*.mat']);
end
Options2.FOVid = cell(1, Ngroups);
NfovsPerGroup = NaN(1, Ngroups);
for k = 1:Ngroups
    NfovsPerGroup(k) = length(Options2.FileLists{1, k});
    Options2.FOVid{k} = cell(1, NfovsPerGroup(k));
    for m = 1:NfovsPerGroup(k)
        Options2.FOVid{k}{m} = Options2.FileLists{1, k}(m).name(end-Options2.lengthFOVid-3:end-4);
    end
end
mkdir([Options2.directory_out, 'Processing2_', Options2.timestamp])
%% Pre-allocate
Results2.FinalTime_min = NaN(max(NfovsPerGroup), Ngroups);
Results2.NfinalROIs = NaN(max(NfovsPerGroup), Ngroups);
Results2.DriftTimeFraction = NaN(max(NfovsPerGroup), Ngroups);
%
Results2.CaT_Frequencies_permin_mean = NaN(max(NfovsPerGroup), Ngroups);
Results2.CaT_Frequencies_permin_median = NaN(max(NfovsPerGroup), Ngroups);
Results2.CaT_Frequencies_permin_cumul = cell(max(NfovsPerGroup), Ngroups);
%
Results2.GDP_Count = NaN(max(NfovsPerGroup), Ngroups);
Results2.GDP_Frequency_permin = NaN(max(NfovsPerGroup), Ngroups);
Results2.GDP_Durations_ms_mean = NaN(max(NfovsPerGroup), Ngroups);
Results2.GDP_Sizes_mean = NaN(max(NfovsPerGroup), Ngroups);
%
%% Extract
fovcount = 0;
for k = 1:Ngroups
    for m = 1:NfovsPerGroup(k)
        %% load and check FOV ID
        fovcount = fovcount + 1;
        fprintf(['FOV# ', num2str(fovcount), ' of ', num2str(sum(NfovsPerGroup)), ' is being processed ... '])
        FileName = [Options2.directory, Options2.subdirectories{1, k}, Options2.FileLists{1, k}(m).name];
        a = load(FileName);
        if strcmp(char(Options2.FOVid{1, k}(m)), a.FOVid); fprintf('FOV ID is valid ... '); else; fprintf('FOV ID is INVALID ... '); end
        %% basics
        Results2.FinalTime_min(m, k) = a.Analysis.FinalTime_min;
        Results2.NfinalROIs(m, k) = a.Analysis.NfinalROIs;
        Results2.DriftTimeFraction(m, k) = sum(a.Events.driftindex)/numel(a.Events.driftindex);
        % *************** CaT frequency ***************
        Results2.CaT_Frequencies_permin_mean(m, k) = mean(cell2mat(a.Analysis.CaT_Frequencies_permin)); % [1/min]
        Results2.CaT_Frequencies_permin_median(m, k) = median(cell2mat(a.Analysis.CaT_Frequencies_permin));  % [1/min]
        if max(cell2mat(a.Analysis.CaT_Frequencies_permin)) > max(Options2.edgesF); error('The maximum CaT frequency exceedes the range defined in "edges"'); end
        aux = histcounts(cell2mat(a.Analysis.CaT_Frequencies_permin), Options2.edgesF, 'Normalization', 'cumcount');
        Results2.CaT_Frequencies_permin_cumul{m, k} = aux ./ max(aux);
        %% GDPs
        Results2.GDP_Count(m, k) = a.Analysis.GDP_Count;
        Results2.GDP_Frequency_permin(m, k) = a.Analysis.GDP_Frequency_permin;
        Results2.GDP_Durations_ms_mean(m, k) = mean(a.Analysis.GDP_Durations_ms); % [ms]
        Results2.GDP_Sizes_mean(m, k) = mean(a.Analysis.GDP_Sizes);
        %% complete
        fprintf('completed\n')
    end
end
%% Figures I: Compute and plot distributions
FieldNames = fieldnames(Results2);
idx = contains(FieldNames, '_cumul');
FieldNames(idx == 0) = [];
for k = 1:numel(FieldNames)
    name = FieldNames{k};
    % get x vector (i.e. starts of bins)
    if contains(name, 'CaT_Frequencies')
        x = Options2.edgesF(1:end-1);
    end
    Distributions2.(name) = NaN(Ngroups*2+1, numel(x));
    Distributions2.(name)(1, :) = x;
    % compute y vectors & plot figure
    figure
    hold on
    LegendNames = cell(1, Ngroups);
    for m = 1:Ngroups
        y = cell2mat(Results2.(name)(:, m));
        n = size(y, 1); % number of rows
        SD = nanstd(y); % SD for each value of x
        SEM = SD ./ sqrt(n); % SEM for each value of x
        MEAN = nanmean(y); % MEAN for each value of x
        Distributions2.(name)(m*2, :) = MEAN;
        Distributions2.(name)(m*2+1, :) = SEM;
        errorbar(x, MEAN, SEM, 'o', ...
            'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black', 'MarkerSize', 1, 'CapSize', 0)
        LegendNames{m} = ['Group ', num2str(m)];
    end
    title(strrep(name, '_', ' ')); legend(LegendNames, 'Location','northeast', 'FontSize', 12)
    ylabel('Cumulative frequency'); xlabel(strrep(name, '_', ' '));
    ylim([0, 1]);
    hold off
    Distributions2.(name) = Distributions2.(name)';
    drawnow;
    saveas(gcf, [Options2.directory_out, 'Processing2_', Options2.timestamp, '\', ...
        'Processing2_', Options2.timestamp, '_', name, '_distr.fig']);
    close(gcf)
end
%% Figures II: mean +/- SEM
FieldNames = fieldnames(Results2);
idx = zeros(numel(FieldNames), 1);
for k = 1:numel(FieldNames)
    name = char(FieldNames(k));
    if strcmp(class(Results2.(name)), 'double')
        idx(k) = 1;
    end
end
idx = logical(idx);
FieldNames = FieldNames(idx);
for k = 1:numel(FieldNames)
    name = char(FieldNames(k));
    plot_meansem(Results2.(name))
    title(strrep(name, '_', ' ')); ylabel(strrep(name, '_', ' ')); xlabel('Measurement groups');
    drawnow;
    saveas(gcf, [Options2.directory_out, 'Processing2_', Options2.timestamp, '\', ...
        'Processing2_', Options2.timestamp, '_', name, '_stats.fig']);
    close(gcf)
end
%% Save MAT files
close all;
clearvars -except Distributions2 NfovsPerGroup Ngroups Options2 Results2
Results2 = orderfields(Results2);
% MAT files
save([Options2.directory_out, 'Processing2_', Options2.timestamp, '\', ...
    'Processing2_', Options2.timestamp, '.mat'])
save([Options2.directory_out, 'Processing2_', Options2.timestamp, '\', ...
    'Processing2_', Options2.timestamp, '_Origin_stats.mat'], '-struct', 'Results2')
save([Options2.directory_out, 'Processing2_', Options2.timestamp, '\', ...
    'Processing2_', Options2.timestamp, '_Origin_distr.mat'], '-struct', 'Distributions2')
%%
disp(['Results were saved to directory ', Options2.directory_out, 'Processing2_', Options2.timestamp, '\'])
warning('on', 'MATLAB:print:FigureTooLargeForPage')
%% Local functions

function plot_meansem(y)
% Plots mean +/- SEM plus individual data points, separately for all
% columns of y
% Input:
% ==========
% y: 2D matrix, where each column is one category
% Output:
% ==========
% figure: graph displaying mean +/- SEM plus individual data points
% ************************************************************
Nrows = size(y, 1);
Ncols = size(y, 2);
n = sum(~isnan(y)); % n per column
SD = nanstd(y); % SD per column
SEM = SD ./ sqrt(n); % SEM per column
MEAN = nanmean(y);
figure
hold on
xpos = 1:Ncols;
er = errorbar(xpos, MEAN, SEM, '-s', 'MarkerSize', 10, 'CapSize', 20, 'LineWidth', 2);
er.Color = [0 0 0]; er.LineStyle = 'none';
xlim([0, size(y, 2)+1]);
xticks(1:Ncols);
xshift = -0.25;
for k = 1:size(y, 2)
    scatter(ones(Nrows, 1)*k + xshift, y(:, k), 'b')
end
hold off

end