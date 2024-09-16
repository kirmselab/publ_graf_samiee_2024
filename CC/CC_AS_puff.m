%% This script computes Membrane Potential (Vm) in CC_AS_puff protocol Recordings
% ==================================================
% Authors: Arash Samiee, Knut Kirmse
% last modification: 07/2024
% ==================================================
% How to use:
% (1) Specify parameters
% (2) Run the script
% (3) Check for Bad Traces
% (4) Transfer the results
% ==================================================
% Inputs:
% (1) Path & File directory
% ==================================================
% Outputs:
% (1) Median Vm per interval [mV]
% ==================================================
%
%% Specify parameters
opt.path = '\folder\';
opt.file = 'file.txt';
opt.fs = 20000; % sampling frequency [Hz]
opt.interval = 10; % [s]
%
%% Clean-up
clearvars -except opt; close all;
%
%% Compute median Vm per sweep
%  read data
data = readmatrix([opt.path, opt.file]);
% some checks
if data(2, 1) ~= 1/opt.fs*1000
    error('Sampling interval is incorrect.')
end
Nintervals = floor(size(data, 1) / (opt.interval * opt.fs));
Nsamples = Nintervals * (opt.interval * opt.fs);
data(Nsamples+1:end, :) = [];
Ihold_max = max(abs(data(:, 3)));
BadTraces_Vm = find(Ihold_max > 1); % 1 pA deviation is accepted
if isempty(BadTraces_Vm)
    disp('No bad traces for computing Median_Vm_perSweep.')
else
    disp(['Bad traces for computing Median_Vm_perSweep: ', num2str(BadTraces_Vm), '.'])
end
% compute quantities
aux = reshape(data(:, 2), Nsamples/Nintervals, []);
Median_Vm_perInterval = median(aux, 1);
clear a b c d Ihold_max aux
%% Save MAT file
%
save([opt.path, opt.file, '_CC_AS_puff.mat'])