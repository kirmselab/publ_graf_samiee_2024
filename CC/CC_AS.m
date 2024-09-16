%% This script computes Membrane Potential (Vm) and Membrane Resistance (Rm) in CC_AS protocol Recordings
% ==================================================
% Authors: Arash Samiee, Knut Kirmse
% last modification: 01/2023
% ==================================================
% How to use:
% (1) Specify parameters
% (2) Run the script
% (3) Check for Bad Traces
% (4) Transfer the results
% ==================================================
% Inputs:
% (1) Path & File directory
% (2) Ihold test Pulse
% ==================================================
% Outputs:
% (1) Rm per each sweep [MOhm]
% (2) Median Vm per each sweep [mV]
% (3) Delta Vm in test pulse period
% (4) Bad Traces for noise control
% ==================================================
%
%% Specify parameters
opt.path = '\folder\';
opt.file = 'file.txt';
opt.fs = 20000; % sampling frequency [Hz]
opt.a = 0; % start of first Vm period (before test pulse) [s]
opt.b = 9.75; % end of first Vm period (before test pulse) [s]
opt.c = 11; % start of second Vm period (after test pulse) [s]
opt.d = 20; % end of second Vm period (after test pulse) [s]
opt.Ihold_testpulse = -20; % Ihold during test pulse [pA]
opt.e = 10; % start of test pulse period for computing Rm [s]
opt.f = 10.25; % end of test pulse period for computing Rm [s]
%
%% Clean-up
clearvars -except opt; close all;
%
%% Compute median Vm per sweep
a = opt.a * opt.fs + 1; % convert time to sample point
b = opt.b * opt.fs; % convert time to sample point
c = opt.c * opt.fs + 1; % convert time to sample point
d = opt.d * opt.fs; % convert time to sample point
%  read data
data = readmatrix([opt.path, opt.file]);
% some checks
if data(2, 1) ~= 1/opt.fs*1000
    error('Sampling interval is incorrect.')
end
Ihold_max = max(abs(data([a:b, c:d], 3:2:end)));
BadTraces_Vm = find(Ihold_max > 1); % 1 pA deviation is accepted
if isempty(BadTraces_Vm)
    disp('No bad traces for computing Median_Vm_perSweep.')
else
    disp(['Bad traces for computing Median_Vm_perSweep: ', num2str(BadTraces_Vm), '.'])
end
% compute quantities
Median_Vm_perSweep = median(data([a:b, c:d], 2:2:end));
clear a b c d Ihold_max
%% Compute membrane resistance per sweep
e = opt.e * opt.fs + 1; % convert time to sample point
f = opt.f * opt.fs; % convert time to sample point
Ihold_max = max(abs(data(e:f, 3:2:end) - opt.Ihold_testpulse));
BadTraces_Rm = find(Ihold_max > 1); % 1 pA deviation is accepted
if isempty(BadTraces_Rm)
    disp('No bad traces for computing Rm_perSweep.')
else
    disp(['Bad traces for computing Rm_perSweep: ', num2str(BadTraces_Rm), '.'])
end
% compute quantities
DeltaVm_perSweep = median(data(e:f, 2:2:end)) - Median_Vm_perSweep;
Rm_perSweep = DeltaVm_perSweep ./ opt.Ihold_testpulse .* 1000;
Median_Vm_perSweep = Median_Vm_perSweep';
DeltaVm_perSweep = DeltaVm_perSweep';
Rm_perSweep = Rm_perSweep';
clear e f Ihold_max data
%% Plot figure and save MAT file
figure
subplot(2, 1, 1)
scatter(1:numel(Median_Vm_perSweep), Median_Vm_perSweep)
ylabel('Median Vm per sweep [mV]')
xlabel('Sweep number []')
subplot(2, 1, 2)
scatter(1:numel(Rm_perSweep), Rm_perSweep)
ylabel('Rm per sweep [MOhm]')
xlabel('Sweep number []')
%
save([opt.path, opt.file, '_CC_AS.mat'])