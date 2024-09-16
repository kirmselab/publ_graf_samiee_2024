%% This script computes 6 parameters for current injection protocol (CC_inj)
% ==================================================
% Authors: Arash Samiee, Knut Kirmse
% last modification: 02/2023
% ==================================================
% How to use:
% (1) Specify parameters
% (2) Run the script
% (3) Visualize the traces and mark any incorrect AP recognition 
% (4) Transfer the results
% ==================================================
% Inputs:
% (1) Path & File directory
% (2) Ihold test Pulse: PW1: (-30:10:120) PW2: (-60:20:240)
% ==================================================
% Outputs:
% (1) Average AP frequency per trace [Hz]
% (2) Median Vm [mV] per Trace before current injection
% (3) AP threshold of the first AP per Trace [mV]
% (4) Maximum instantaneous AP frequency [Hz] per trace
% (5) Rebound depolarization in responce to current change [mV]
% (6) Voltage Sag in responce to negative current injection [mV]
% ==================================================
%
%% Specify parameters
opt.path = '\folder\';
opt.file = 'file.txt';
opt.fs = 20000; % sampling frequency [Hz]
opt.t1 = 1; % start of current injection [s]
opt.t2 = 1.5; % end of current injection [s]
opt.Ihold_testpulse = -30:10:120; % delta Ihold during current injection [pA]
opt.smoothWindow = 5; % window length [sample points] of Savitzky-Golay smoothing (applied to dVm/dt)
opt.smoothOrder = 2; % order [] of Savitzky-Golay smoothing  (applied to dVm/dt)
opt.dVm_dt_Threshold = 20000; % dVm/dt threshold for AP detection [mV/s]
opt.minAPdelay = 2; % APs with an Inter-AP-Interval lower than this value [ms] will be deleted
opt.smoothWindowSag = 400; % window length [sample points] of Savitzky-Golay smoothing (voltage sag measurment ONLY)
opt.searchWindowSag = 0.1; % window length [s] for searching min/max of voltage sage and rebound depolarization
%
%% Clean-up
clearvars -except opt; close all;
%
%% Compute
Ntraces = numel(opt.Ihold_testpulse);
t1 = opt.t1 * opt.fs + 1; % convert time to sample point
t2 = opt.t2 * opt.fs; % convert time to sample point
%  read data
data = readmatrix([opt.path, opt.file]);
Nsamples = size(data, 1);
% some checks
if data(2, 1) ~= 1/opt.fs*1000
    error('Sampling interval is incorrect.')
end
if Ntraces ~= (size(data, 2)-1) / 2
    error('Number of traces in data does not fit.')
end
% compute quantities
dVm_dt = NaN(Nsamples, Ntraces); % pre-allocate
AboveThreshold = false(Nsamples, Ntraces);
APonsets = cell(1, Ntraces); % AP onsets within interval [t1; t2][s]
APthresholds = cell(1, Ntraces); % AP thresholds [mV]
FirstAPthresholds = NaN(Ntraces, 1); % AP threshold of the first AP per Trace [mV]
NumberAPperTrace = NaN(1, Ntraces);
MaxFrequperTrace = NaN(Ntraces, 1); % maximum instantaneous AP frequency [Hz] per trace
BaselineVm = NaN(Ntraces, 1); % median Vm [mV] per Trace before current injection
for k = 1:Ntraces
    aux = gradient(data(:, k*2), 1/opt.fs); % 1st derivative
    dVm_dt(:, k) = smooth(aux, opt.smoothWindow, 'sgolay', opt.smoothOrder);
    AboveThreshold(:, k) = dVm_dt(:, k) > opt.dVm_dt_Threshold;
    aponsets = find(diff(AboveThreshold(:, k)) == 1) + 1; % AP onset of k-th trace [samples]
    aponsets = aponsets(aponsets > t1 & aponsets <= t2);
    if numel(aponsets) > 1
        delays = find(diff(aponsets) < opt.fs/1000*opt.minAPdelay) + 1;
        aponsets(delays) = [];
    end
    NumberAPperTrace(1, k) = numel(aponsets)';
    APonsets{1, k} = (aponsets-1) ./ opt.fs; % AP onsets [s]
    APthresholds{1, k} = data(aponsets, k*2);
    if NumberAPperTrace(1, k) >= 1
        FirstAPthresholds(k, 1) = data(aponsets(1), k*2);
    end
    if NumberAPperTrace(1, k) >= 2
        MaxFrequperTrace(k, 1) = 1/min(diff(APonsets{1, k}));
    end
    BaselineVm(k, 1) = median(data(1:t1-1, k*2));
end
AvgFrequperTrace = NumberAPperTrace ./ (opt.t2 - opt.t1); % mean AP frequency [Hz] per trace
AvgFrequperTrace = AvgFrequperTrace';
%
negSteps = find(opt.Ihold_testpulse < 0); % steps with negative current injection
VoltageSag = NaN(numel(negSteps), 1);
ReboundDepol = NaN(numel(negSteps), 1);
for k = negSteps
    aux1 = data(:, k*2);
    aux2 = smooth(aux1, opt.smoothWindowSag, 'sgolay', opt.smoothOrder);
    VoltageSag(k, 1) = min(aux2(t1:t1+opt.searchWindowSag*opt.fs-1)) - median(aux2(t2-opt.searchWindowSag*opt.fs:t2));
    ReboundDepol(k, 1) = max(aux2(t2+1:t2+opt.searchWindowSag*opt.fs)) - BaselineVm(k, 1);
end
%% Plot in Viewer app
CC_INJ_TracesViewer

%%  Save MAT file
clear t1 t2 aponsets k aux aux1 aux2 AboveThreshold data delays dVm_dt negSteps
%save([opt.path, opt.file, '_CC_inj.mat'])