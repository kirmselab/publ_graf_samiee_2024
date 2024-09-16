%% This script runs UFARSA_OGB1
% ----------------------------------------------------------------------
% Purpose:
% The algorithm applies UFARSA (Rahmati et al. J Neurophysiol (2018)).
% ======================================================================
% Use:
% (1) Set parameters in section "Specify parameters".
% (2) Optional: Set further parameters of UFARSA
%     (stored in UFARSA_master\internal_parameters.m.).
% (3) >> runUFARSA_OGB1
% ======================================================================
% Help:
% N/A.
% ======================================================================
%
%% Specify parameters
clear variables; close all;
% ******************************************************************************************************
% Parameters I: General
% ******************************************************************************************************
options.firstframe = 0;         % Set to 1, if the first frame in the TIF stack is NOT the first frame acquired,
%                                 but e.g. an average used for stack alignment. Otherwise set to 0. (by default 0)
options.smooth = 0;             % 1: Apply Savitzky-Golay smoothing to ROI-wise and pixel-wise F(t). 0 - Do not apply. (by default 0)
options.SG_length = 18;         % Length [frames] of Savitzky-Golay smoothing window applied to (1) pixel-wise F(t) and (2) 1st derivative of F_perROI
options.SG_order = 2;           % Order [] of Savitzky-Golay smoothing
options.window_F0 = 2000;       % Length [frames] of Savitzky-Golay smoothing window for computing F0(t). [] – use all frames.
options.F0_correct = -0.4;      % If the minimum DF/F0(t) for a given ROI is lower than options.F0_correct,
%                                 F0(t) is instead calculated as the percentile specified by options.F0_percentile.
options.F0_percentile = 20;     % If the minimum DF/F0(t) is lower than options.F0_correct,
%                                 F0(t) is instead calculated as the percentile specified by options.F0_percentile.
options.ROIexpansion = 0;       % Euclidean distance [pixels] by which ROIs will be expanded. 0 – no expansion.
options.SaveOptVars = 1;        % 1: Save optional variables. 0: Skip this step. (by default 0)
options.substractOffset = 0;    % This value will be subtracted vom all data points in F(t).
% ******************************************************************************************************
% Parameters II: UFARSA (for the rest of parameters, see "internal_parameters.m" file)
% ******************************************************************************************************
opt.scale_NoiseSTD = 2;         % Leading-threshold scaling constant (by default 2.5)
opt.remove_drifts = 1;          % 1: Remove slowly varying drifts. 0: Skip the drift removal step. (by default 1)
opt.remove_posDeflections = 0;  % 1: Apply large-impulse (deflection) removal step. 0: Skip this step. (by default 0)
opt.remove_negDeflections = 1;  % 1: Remove large short-lasting negative deflections. 0: Skip this step. (by default 1)
opt.demerging = 1;              % 1: Apply the demerging step. 0: Skip this step. (by default 1)
%
%% Get TIF stack
options.timestamp = datestr(now, 'yymmdd_HHMMSS');
[options.FileNameTIF,options.PathNameTIF] = uigetfile('*.tif', 'Select the TIF file');
if isequal(options.FileNameTIF,0)
   disp('User selected Cancel')
   return;
else
    [~,options.FileNameOut,~] = fileparts([options.PathNameTIF, options.FileNameTIF]);
    disp(['TIF stack: ', fullfile(options.PathNameTIF, options.FileNameTIF)])
end
% aux = strfind(options.PathNameTIF, '\');
% prefix = options.PathNameTIF(aux(end-2)+1:end);
% prefix = strrep(prefix, '\', '_');
% options.FileNameOut = [prefix, options.FileNameOut];
data = ReadTIFFobj([options.PathNameTIF, options.FileNameTIF]);
data = double(data);
data = data - options.substractOffset;
if options.firstframe == 1
    data(:, :, 1) = [];
end
Nrows = size(data, 1);
Ncols = size(data, 2);
Nframes = size(data, 3);
%
%% Get TXT with ROIs and optionally expand ROIs
[options.FileNameROIs,options.PathNameROIs] = uigetfile('*.txt', 'Select the TXT file that specifies the ROIs', options.PathNameTIF);
if isequal(options.FileNameROIs,0)
   disp('User selected Cancel')
   return;
else
   disp(['ROI list TXT file: ', fullfile(options.PathNameROIs, options.FileNameROIs)])
end
coord = dlmread([options.PathNameROIs, options.FileNameROIs], '\t', 1, 0);
Nrois = 0.5 * size(coord, 2);
% Store linear indices of initial ROIs
coord_ind = cell(1, Nrois);
for k = 1:Nrois
    ind = sub2ind([Nrows, Ncols], coord(:, k*2-1), coord(:, k*2));
    ind(isnan(ind)) = [];
    coord_ind{1, k} = ind;
end 
% Expand ROIs (if applicable)
if options.ROIexpansion > 0
    coord_extended = NaN(Nrows * Ncols, size(coord, 2));
    for k = 1:Nrois
        im = false(Nrows, Ncols);
        ind = coord_ind{1, k};
        im(ind) = 1;
        im = bwdist(im) <= options.ROIexpansion; % Euclidean distance transform of the binary ROI image
        ind = find(im == 1);
        [rows,cols] = ind2sub([Nrows, Ncols],ind);
        coord_extended(1:length(rows), k*2-1) = rows;
        coord_extended(1:length(cols), k*2) = cols;
    end
    coord_extended(all(isnan(coord_extended), 2), :) = [];
else
    coord_extended = coord;
end
% Store linear indices of extended ROIs
coord_extended_ind = cell(1, Nrois);
for k = 1:Nrois
    ind = sub2ind([Nrows, Ncols], coord_extended(:, k*2-1), coord_extended(:, k*2));
    ind(isnan(ind)) = [];
    coord_extended_ind{1, k} = ind;
end 
clear k im rows cols ind
%% Get TXT with driftperiods
[options.FileNameDrift,options.PathNameDrift] = uigetfile('*.txt', 'Select the TXT file that specifies the drift periods', options.PathNameTIF);
if isequal(options.FileNameDrift,0)
   disp('User selected Cancel')
   return;
else
   disp(['Drift periods TXT file: ', fullfile(options.PathNameDrift, options.FileNameDrift)])
end
s = dir([options.PathNameDrift, options.FileNameDrift]);
if s.bytes == 0
    driftindex = false(1, Nframes);
else
    driftperiods = dlmread([options.PathNameDrift, options.FileNameDrift], '\t', 0, 0);
    % Convert 'driftperiods' to logical 'driftindex'
    driftindex = false(1, Nframes);
    for k = 1:size(driftperiods, 1)
        for j = driftperiods(k, 1):driftperiods(k, 2)
            driftindex(1, j) = 1;
        end
    end
end
clear j s k
% Compute final frames for template matching
finalFrames = find(driftindex == 0);
%
%% Get TXT with file onset frames
[options.FileNameFileOnsetFrames,options.PathNameFileOnsetFrames] = uigetfile('*.txt', 'Select the TXT file that specifies the file onset frames', options.PathNameTIF);
if isequal(options.FileNameFileOnsetFrames,0)
   disp('User selected Cancel')
   return;
else
   disp(['File onset frame list TXT file: ', fullfile(options.PathNameFileOnsetFrames, options.FileNameFileOnsetFrames)])
end
FileOnsetFrames = dlmread([options.PathNameFileOnsetFrames, options.FileNameFileOnsetFrames], '\t', 0, 0);
%
%% Compute F(t), F0(t) and DF/F0(t) per ROI
% ROI expansion is applied here.
disp('Computing ...')
F_perROI = NaN(Nrois, Nframes);
F0_perROI = NaN(Nrois, Nframes);
DF_F0_perROI = NaN(Nrois, Nframes);
Npixels = Nrows * Ncols;
for k = 1:Nrois
    ind = coord_extended_ind{1, k};
    for m = 1:Nframes
        F_perROI(k, m) = nanmean(data(ind + ((m-1) * Npixels)));
    end
    if options.smooth == 1
        F_perROI(k, :) = smooth(F_perROI(k, :), options.SG_length, 'sgolay', options.SG_order);
    end
    if isempty(options.window_F0) == 1
        F0_perROI(k, :) = median(F_perROI(k, :), 'omitnan');
    else
        F0_perROI(k, :) = movmedian(F_perROI(k, :), options.window_F0, 'omitnan');
    end
    if sum(driftindex) > 0
        F0_perROI(k, driftindex) = NaN;
    end     
    DF_F0_perROI(k, :) = bsxfun(@rdivide, F_perROI(k, :), F0_perROI(k, :)) - 1;
end
% Correct F0(t) calculation for 'highly active' cells
DF_F0_perROI_min = min(DF_F0_perROI, [], 2);
ROIs_Prctile_F0 = zeros(1, Nrois);
for k = 1:Nrois
    if DF_F0_perROI_min(k) < options.F0_correct
        ROIs_Prctile_F0(1, k) = 1;
        F0_perROI(k, :) = prctile(F_perROI(k, :), options.F0_percentile, 2);
        if sum(driftindex) > 0
            F0_perROI(k, driftindex) = NaN;
        end
        DF_F0_perROI(k, :) = bsxfun(@rdivide, F_perROI(k, :), F0_perROI(k, :)) - 1;
    end
end
%
DF_F0_perROI = DF_F0_perROI';
clear ind k m Npixels Nrows Ncols coord_ind coord_extended_ind F0_perROI
%
clear k
%
%% Prepare cell array 'values'
values = cell(6, Nrois);
for k = 1:Nrois
    values{1, k} = ['ROI_', int2str(k)];
end
% In CATHARSiS, values(2:4, :) contain F(xyt), F0(xyt) and DF(xyt) and
% hence remain unused here.
% In CATHARSiS, values(5, :) contain candidate DF(xy) templates and hence remain
% unused here.
% In CATHARSiS, values(6, :) contain final DF(xy) templates and hence remain
% unused here.
clear data
clear k
finalROIs = 1:Nrois;
non_finalROIs = setdiff(1:Nrois, finalROIs);
% In CATHARSiS, ROIs without final template will be excluded.
%
%% Compute DETECTION CRITERION
% For clarity, we here keep the code structure used in CATHARSiS despite
% the redundancy.
DETECTION_CRITERION = NaN(Nframes, Nrois);
for k = finalROIs
    DETECTION_CRITERION(:, k) = DF_F0_perROI(:, k);
end
clear k
%
%% Apply UFARSA to DETECTION CRITERION
opt.which_ROIs = [];
opt.gen_FR_count = 0;       % 1: generate the estimated firing rate vector based on the reconstructed spike-count train, 0: skip it (by default 0)
opt.gen_FR_count_dem = 0;   % 1: generate the estimated firing rate vector based on the reconstructed demerged spike-count train, 0: skip it (by default 0)
opt.FluorFile_name = 'none.non'; % placeholder for UFARSA (see run_UFARSA.m: *_*)
opt.FluorFile_dir  = options.PathNameTIF;
opt.samples_randomly = 0; % parameter in UFARSA: samples for estimating noise STD set in a non-random manner
% parameter in UFARSA: opt.nSamplesForSTD was re-defined, see also changes in 'smoothing_UFARSA.m'
opt.nSamplesForSTD = find(driftindex == 0);
UFARSA_output = cell(1, Nrois);
DETECTION_CRITERION(isnan(DETECTION_CRITERION)) = 0; % UFARSA does not accept NaNs
for k = finalROIs
    [output_UFARSA, opt_out, ~] = run_UFARSA(opt, DETECTION_CRITERION(:, k));
    UFARSA_output{k}.output_UFARSA = output_UFARSA;
    UFARSA_output{k}.opt_out = opt_out;
end
clear k output_UFARSA opt_out
%
%% Correct drift periods
% UFARSA can generate false positive event onsets at first frames of files and/or
% non-drift periods if there is a rise/decay of a CaT whose actual onset is before that frame.
% We therefore set all these frames to 1 in 'driftindex' and re-define 'finalFrames' accordingly.
%
% add first frame after each drift period to driftindex
driftindex = movmax(driftindex, [1, 0]);
% add file onset frames to driftindex
driftindex(FileOnsetFrames') = 1;
finalFrames = find(driftindex == 0);
%
%% Output TXT (to be displayed with Clampfit)
Npara = 5; % number of output parameters
% (I)   UFASAR-demerged event times
% (II)  UFASAR-smoothed DETECTION CRITERION
% (III) DETECTION CRITERION
% (IV)  DF/F0
% (V)   drift periods
Output = zeros(Nframes, Npara * Nrois); % Clampfit cannot handle NaNs.
for k = finalROIs
    Output(:, (k-1) * Npara + 1) = UFARSA_output{k}.output_UFARSA.eTrain_dem';
    Output(:, (k-1) * Npara + 2) = UFARSA_output{k}.output_UFARSA.fluors.afterSmoothing';
    Output(:, (k-1) * Npara + 3) = UFARSA_output{k}.output_UFARSA.fluors.original';
    Output(:, (k-1) * Npara + 4) = DF_F0_perROI(:, k);
end
%
if sum(driftindex) > 0
    % Set output signals I to IV to zero for all drift periods
    Output(driftindex', :) = 0;
    % Add output signal V) drift periods
    for k = finalROIs
        Output(:, (k-1) * Npara + 5) = driftindex';
    end
end
dlmwrite([options.PathNameTIF, options.FileNameOut, '_UFARSAevents_DCsmoothed_DC_DFF0_Drift_', options.timestamp, '.txt'], Output, '\t');
disp(['Traces to be displayed with Clampfit were exported to: ', options.PathNameTIF, options.FileNameOut, '_UFARSAevents_DCsmoothed_DC_DFF0_Drift_', options.timestamp, '.txt']);
%
%% Save MAT file
% Group results into table
rownames = {'ROI number', 'F(xyt)', 'F0(xyt)', 'DF(xyt)', 'Rise peaks', 'Candidate templates (xyn)', 'Unused', 'User-selected templates (0 - deselected)', 'Final template (xy)', 'DF/F0(t) per ROI w/o expansion', 'Detection criterion', 'UFARSA: Demerged event train', 'UFARSA: leading_thr', 'UFARSA: smoothing_param'};
Results_perROI = cell2table(cell(14, Nrois), 'rownames', rownames);
Results_perROI{[1:4, 6, 9], 1:Nrois} = values(1:6, :);
Results_perROI{10, finalROIs} = num2cell(DF_F0_perROI(:, finalROIs), 1);
Results_perROI{11, finalROIs} = num2cell(DETECTION_CRITERION(:, finalROIs), 1);
clear rownames values DF_F0_perROI DETECTION_CRITERION
% Note: Only values of 'output_UFARSA.eTrain_dem'' in row 12 of
% 'Results_perROI' are set to NaN according to final 'driftindex'.
for k = finalROIs
    aux = UFARSA_output{1, k}.output_UFARSA.eTrain_dem';
    % Set values within drift periods back to NaN
    % Set values at file frame onsets to NaN as UFARSA can generate false
    % positive events
    aux(driftindex') = NaN;
    Results_perROI{12, k} = num2cell(aux, 1);
    Results_perROI{13, k} = num2cell(UFARSA_output{1, k}.output_UFARSA.leading_thr);
    Results_perROI{14, k} = num2cell(UFARSA_output{1, k}.opt_out.smoothing_param);
end
% Save
clear k aux Npara
if options.SaveOptVars == 0
    clear coord coord_extended DF_F0_perROI_min F_perROI ROIs_Prctile_F0
end
save([options.PathNameTIF, options.FileNameOut, '_', options.timestamp, '_events.mat']);
disp(['Matlab workspace was saved to: ', options.PathNameTIF, options.FileNameOut, '_', options.timestamp, '_events.mat']);
disp('Finished.')
%