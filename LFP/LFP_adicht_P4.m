%% LFP_adicht.m
% Author: Jürgen Graf
% date: 21.07.2021
% intended use: extract LFP and Respiration/Movement data from LabChart files (.adicht)
% specifically written for LFP recordings in Gi-DREADD mice
% Prerequisite: 
%           - import_adicht.m (J. Graf)
%           - Simple Data File SDK (from ADInstruments)
%           - Chronux toolbox
%           - peakfinder (Nathanael C. Yoder)
% features: - import data directly from .adicht files
%           - specify recording sections and LFP bands for analysis
%           - downsample to desired sampling rate
%           - calculate Power spectrum using pwelch, also for subsections of desired length
%           - PSD over time using mtspecgram (Chronux toolbox)
%           - event detection parameters is optimized for P3-4 (network
%           bursts showing a power increase in 15-30 Hz (beta) range)

%% Get data
% extract metadata
clear; close all
format long
primaryFolder = 'D:\JG\Project Gi-DREADD';
[options.file options.file_path] = uigetfile({'*.*',  'All Files (*.*)'}, 'Select a file', primaryFolder);
[a, options.filename, options.file_ext] = fileparts([options.file_path options.file]); % a is again the file_path
clear a
options.timestamp = datestr(now, 'yymmdd_HHMMSS');
[f, data, comments] = import_adicht([options.file_path options.file]);
clear primaryFolder
%% Print comments
comments
%% specify time windows for control and C21 recording
options.LFP_Ref =  [1729 2132]; % Reference recording for LFP, when electrode was still outside the brain
options.control = [3675 5613]; % control recording in s from start of file
% how to deal with different records?
options.C21 = [5690 15543]; % treatment recording in s from start of file
cond = {'control', 'C21', 'LFP_Ref'};
cond_color = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250]; % specific colors for each condition
for n = 1:length(cond)
    field = string(cond{n});
    bands.(field){1,1} = [4 8]; % theta band
    bands.(field){2,1} = [8 12]; % alpha band
    bands.(field){3,1} = [15 30]; % beta band
    bands.(field){4,1} = [4 15]; % custom data driven frequency band (P10)
    bands.(field){5,1} = [8 25]; % custom data driven frequency band for HBOs (P4)
    bands.(field){6,1} = [10 16]; % custom data driven frequency band for HBOs (P4), even more specific    
end
drift_ext = 0.10; % extension in seconds to remove after noise peak detection; should not be smaller than 1/lfp.fs
nSec_base = 30; % window length in seconds for baseline calculation
SPW_thr = -0.1; % LFP threshold for SPW detection
std_factor = 12; % factor for determination of threshold (median + factor * MAD)
nSec_win = 1; % length of sliding window for spectrogram calculation in sec
win_overlap = 0.8; % overlap in %
win = [nSec_win nSec_win-win_overlap*nSec_win];  % length and overlap of sliding window
Spec_fs = 1/(nSec_win-nSec_win*win_overlap); % temporal resolution of spectrum in Hz, depending on window size and overlap
MinEventTime = 0.6; % Minimum time an event is showing a power over threshold
metadata.resp.Ch = find(contains(f.channel_names,'Move')); % find Respiration/Movement channel
metadata.lfp.Ch = find(strncmpi(f.channel_names,'Ephys', 5)); % find Ephys channel
metadata.lfp.origFs = f.channel_specs(1, metadata.lfp.Ch).fs(1); % original recording frequency of LFP in Hz
metadata.resp.origFs = f.channel_specs(1, metadata.resp.Ch).fs(1); % original recording frequency of Resp in Hz
metadata.lfp.fs = 500; % specify target lfp fs
metadata.lfp.fsRatio = metadata.lfp.origFs/metadata.lfp.fs; % for resample this must be an integer
metadata.resp.fs = 100; % target resp fs
metadata.resp.fsRatio = metadata.resp.origFs/metadata.resp.fs;
metadata.nRecords = f.n_records;
metadata.FileStartTime = f.records.data_start; % Datenumber
metadata.TimeofDay =  datestr(metadata.FileStartTime,'HH:MM:SS PM');
metadata.comments = comments; % list of comments in LabChart
if range(f.channel_specs(1, metadata.lfp.Ch).fs) ~= 0
    error('Recording frequency for lfp was changed during recording!')
end
if range(f.channel_specs(1, metadata.resp.Ch).fs) ~= 0
    error('Recording frequency for Respiration was changed during recording!')
end
%% open previous .mat file to import time windows
[options.OldFileName,options.OldPathName] = uigetfile('*.mat', 'Select .mat file to import time windows from', options.file_path);
S = load([options.OldPathName, options.OldFileName], 'options');
options.LFP_Ref =  S.options.LFP_Ref; % Reference recording for LFP, when electrode was still outside the brain
options.control = S.options.control; % control recording in s from start of file
options.C21 = S.options.C21; % treatment recording in s from start of file
clear S
%% truncate data in dependence of selected time windows for control and C21
for n = 1:length(cond)
    field = string(cond{n});
    bands.(field)(:,2:end) = [];
end
clear driftidx driftperiods driftperiods_all eventperiods goodidx lfp mt_s mt_s_log mt_t osci pxx SPW subsec w z_mt_s
for n = 1:size(cond,2) % control, C21 and LFP reference
    field = string(cond{n});
    lfp{1,n} = data{1,metadata.lfp.Ch}(options.(field)(1)*metadata.lfp.origFs:options.(field)(2)*metadata.lfp.origFs);
    metadata.rectime{1, n} = length(lfp{1,n})/metadata.lfp.origFs; % total recording time in s
end
if ~isempty(data{:,metadata.resp.Ch})
    resp{1,1} = data{1,metadata.resp.Ch}(options.control(1)*metadata.resp.origFs:options.control(2)*metadata.resp.origFs);%Respiration/Movement signal under control
    resp{1,2} = data{1,metadata.resp.Ch}(options.C21(1)*metadata.resp.origFs:options.C21(2)*metadata.resp.origFs);%Respiration/Movement signal under C21
end

%% Get TXT with driftperiods
[options.FileNameDrift,options.PathNameDrift] = uigetfile('*.txt', 'Select the TXT file that specifies drift periods', options.file_path);
% tab delimited table with drift start (column 1) and end (column 2) in seconds
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
end
clear s
%%  Detect lfp noise and delete drift periods 
tic
for n = 1:size(lfp,2) % control, C21 and LFP reference
    field = string(cond{n});
    lfp{3,n} = single([diff(lfp{1,n});0]); % First derivative of raw LFP
    driftidx{1,n} = [find(lfp{3,n}>0.17); find(lfp{3,n}<-0.17); find(abs(diff(lfp{3,n}))>0.3)]; % detect high amplitude noise peaks, meaning 0.2 mV per sample point, which is a really sharp increase
    lfp{3,n} = []; % delete first derivative of lfp for saving space
    lfp{4,n} = false(length(lfp{1,n}), 1); % create logical vector 
    driftidx_ext{1,n} = unique(reshape(cumsum([driftidx{1,n},ones(length(driftidx{1,n}),metadata.lfp.origFs*drift_ext)],2),1,[]), 'sorted'); % generate driftindex (drift = 1) and extend all detected drifts by drift_ext (sec)
    lfp{4,n}(driftidx_ext{1,n}) = 1;
    for k = 1:size(driftperiods, 1)
        driftStart = int64((driftperiods(k, 1)-options.(field)(1))*metadata.lfp.origFs);
        driftEnd = int64((driftperiods(k, 2)-options.(field)(1))*metadata.lfp.origFs);
        if driftStart < 1 || driftStart > length(lfp{4,n})
            continue
        else
            lfp{4,n}(driftStart:driftEnd) = 1;
        end
    end
    driftperiods_all{1,n}(:,1) = find(diff(lfp{4,n}) == 1);
    driftperiods_all{1,n}(:,2) = find(diff(lfp{4,n}) == -1);
    driftperiods_all{1,n} = driftperiods_all{1,n} ./ metadata.lfp.origFs + options.(field)(1);
    lfp{1,n}(lfp{4,n} == 1) = []; % remove drift from lfp
    driftidx{2,n} = single(find(lfp{4,n} == 1)); % drift index after drift extension and txt with driftperiods
    goodidx{1,n} = single(find(lfp{4,n} == 0)); % index of lfp values outside drift, used to reconstruct original event times
    lfp{4,n} = [];
end
clear n k driftidx_ext driftStart driftEnd
toc
%% Resample and cast to single
% for respiration and movement signal
if ~isempty(data{:,metadata.resp.Ch}) %only resample when resp was imported
    for n = 1:length(resp) % control and C21
        resp{1,n} = single(resample(resp{1,n},1,metadata.resp.fsRatio)); % resampling at 1/fsRatio times original sampling rate, cast to single
    end
end
% for lfp
for n = 1:size(lfp,2) % control, C21 and LFP reference
    lfp{1,n} = single(resample(lfp{1,n},1,metadata.lfp.fsRatio)); % resample lfp
%     lfp{4,n} = lfp{4,n}(1:metadata.lfp.fsRatio:end); % downsample lfp noise vector the standard way
%     lfp{4,n} = movmax(lfp{4,n}, 10); % extend noise peaks again, why?
%     driftidx{2,n} = find(lfp{4,n} == 1); % drift index after resample
%     goodidx{1,n} = single(find(lfp{4,n} == 0)); % index of lfp values outside drift, used to reconstruct original event times
%     lfp{4,n} = []; % delete lfp noise vector
%     lfp{1,n}(driftidx{2,n}) = []; % remove drift from resampled lfp
    metadata.rectime{2, n} = length(lfp{1,n})/metadata.lfp.fs; % total recording time in s, without drift
    lfp{2,n}(:,1) = single([1/metadata.lfp.fs:1/metadata.lfp.fs:metadata.rectime{2,n}]); % time vector, without drift
    lfp{3,n} = lfp{1,n}-medfilt1(lfp{1,n},metadata.lfp.fs); % Baseline correction with median filter width of 1 s
end
% clear data
%% PSD of entire recording using Welch's method
win_welch = 6; % window length in sec for Welch's PSD estimate
overlap = 0.2; % window overlap, number between 0 and 1
win_length = win_welch*metadata.lfp.fs; % larger windows increase resolution for smaller frequency ranges
noverlap = round(overlap*win_length); % number of data points for window overlap
% for lfp
for n = 1:size(lfp,2) % control, C21 and LFP reference
    [pxx{n}, w] = pwelch(lfp{1,n}, win_length, noverlap,[], metadata.lfp.fs, 'psd');   
end
% calculate Bandpower for lfp
for n = 1:size(lfp,2) % control, C21 and LFP reference
    field = string(cond{n});
    for i = 1:size(bands.(field),1)
        bp.(field)(i,1) = bandpower(pxx{n}, w, bands.(field){i,1} , 'psd'); % Bandpower for all specified bands
    end
end

figure
hold on
for n = 1:size(lfp,2) % control, C21 and LFP reference
    plot(w, pxx{n}, 'color', cond_color(n,:))
end
title('PSD using pwelch')
legend(cond, 'Interpreter', 'none')
xlabel('Frequency [Hz]')
ylabel('Power')
xlim([3 40])
hold off
export_fig(([options.file_path, options.filename, '_', options.timestamp]), '-pdf', '-silent', '-transparent', '-nocrop')
%% PSD of subsections of the recording to see if power increases monotonically during recording
subsec.length = 300; % length of subsections in seconds
for n = 1:size(lfp,2) % control, C21 and LFP reference
    subsec.n{n} = floor(metadata.rectime{2, n}/subsec.length); % number of possible subsections using the specified length
    field = string(cond{n});
    for i = 1:subsec.n{n} % for the number of subsections
        o = (i-1)*subsec.length*metadata.lfp.fs + 1; % beginning of subsection
        p = i*subsec.length*metadata.lfp.fs; % end of subsection
        subsec.lfp.(field){i} = lfp{1,n}(o:p); % extract subsection
        subsec.pxx.(field){i} = pwelch(subsec.lfp.(field){i}, win_length, noverlap,[], metadata.lfp.fs, 'psd'); % calculate PSD of subsection
        subsec.bp.(field){i} = bandpower(subsec.pxx.(field){i}, w, bands.control{5, 1} , 'psd'); % calculate bandpower of subsection, bands.control{5, 1} = 8-25 Hz
    end
end
clear o p 
%% plotting PSD and Bandpower for subsections
cmap = colormap('autumn');
figure('Position',[300 100 900 1000])
tl = tiledlayout('flow','TileSpacing','compact'); 
for n = 1:size(lfp,2) % control, C21 and LFP reference
    field = string(cond{n});
    c = (1:floor(255/(subsec.n{n}-1)):256);
    nexttile
    hold on
    for i = 1:subsec.n{n} % for the number of subsections
        plot_color = cmap(c(i), :);
        plot(w, subsec.pxx.(field){i}, 'Color', plot_color)
    end
    title(cond{n}, 'Interpreter', 'none');
    xlabel('Frequency [Hz]')
    ylabel('Power')
    xlim([3 40])
    legend('off')
    hold off

    nexttile
    if subsec.n{n} ~= 0
        b = bar(1:subsec.n{n}, cell2mat(subsec.bp.(field)), 'FaceColor', 'flat');
        for i = 1:subsec.n{n} % for the number of subsections
            plot_color = cmap(c(i), :);
            b.CData(i, :) = plot_color;
        end
    end
    title(cond{n}, 'Interpreter', 'none');
    xlabel('Subsections')
    ylabel('Bandpower [8-25 Hz]')
end
title(tl, ['PSD of subsections of length ', num2str(subsec.length), ' s'])
export_fig(([options.file_path, options.filename, '_', options.timestamp]), '-pdf', '-silent', '-transparent', '-nocrop', '-append') %-preserve_size
%% PSD using mtspecgramc from chronux toolbox
params.fpass = [3 100]; % frequency band to be used in the calculation
params.tapers = [1 2]; % [TW K] TW is time-bandwidth product, K is number of tapers, K<=(TW*2)-1
params.trialave = 0; % 1 - average over trials, 0 - don't average
params.err = [1 .05]; % error calculation, [1 p] - Theoretical error bars; [2 p] - Jackknife error bars; [0 p] or 0 - no error bars) - optional. Default 0.
params.Fs = metadata.lfp.fs; % sampling frequency
params.pad = 1;% if number of samples n = 500, pad = -1, no padding, pad = 0, we pad the fft to 512 bins (take to the next nearest power of 2^x)
% if pad = 1, we pad to 1024 bins, if pad = 2, we pad to 2048 bins
for n = 1:size(lfp,2) % control, C21 and LFP reference
    [mt_s{n}, mt_t{n}, mt_f, serr{n}] = mtspecgramc(lfp{1,n},win,params); % 
    mt_s{n} = mt_s{n}'; % Power/Frequency, dB/Hz
    mt_s_log{n} = 20*log10(abs(mt_s{n})); % conversion to dB as RMS-spectrum (root-mean-square spectrum)
    z_mt_s{n} = zscore(mt_s{n},0,2);
    mt_t{n} = mt_t{n} - 0.102; % correction to set center of first bin to 0.5 sec
end
clear serr
%% Calculating bandpower and a dynamic threshold depending on a moving window for baseline calculation
for n = 1:size(lfp,2) % control, C21 and LFP reference
    field = string(cond{n});
    for i = 1:size(bands.(field),1)
        bands.(field){i,2}(:,1) = bandpower(mt_s{n}, mt_f, bands.(field){i,1}, 'psd'); % absolute bandpower across time
        bands.(field){i,3} = movmedian(bands.(field){i,2}, nSec_base*Spec_fs); % baseline power as median over a window of length nSec_base
        bands.(field){i,4} = movmad(bands.(field){i,2}, nSec_base*Spec_fs*3); % MAD for a time window three times larger
        bands.(field){i,5} = bands.(field){i,3}+bands.(field){i,4}.*std_factor; % threshold for detection of discrete events
    end
end
%% Plotting the entire recording for visual inspection
prompt = 'What condition do you want to plot? 1:control; 2:C21: 3:LFP_Ref: ';
n = input(prompt);
field = string(cond{n});
fprintf('1: theta %d-%d Hz\n',bands.control{1, 1})
fprintf('2: alpha %d-%d Hz\n',bands.control{2, 1})
fprintf('3: beta %d-%d Hz\n',bands.control{3, 1})
fprintf('4: theta-alpha %d-%d Hz\n',bands.control{4, 1})
fprintf('5: custom HBO %d-%d Hz\n',bands.control{5, 1})
fprintf('6: custom HBO %d-%d Hz\n',bands.control{6, 1})
prompt = 'What bandpower do you want to plot?: ';
i = input(prompt);
bandlabel = ['[', num2str(bands.control{i, 1}(1)),'-',num2str(bands.control{i, 1}(2)), ' Hz]'];
freqs = round(logspace(log10(1),log10(100),100),2);
figure('Position',[250 50 1200 700]) % [left bottom width height]
tl = tiledlayout(3,1,'TileSpacing','compact','Padding','none'); 
ax1 = nexttile;
plot(lfp{2,n},lfp{1,n});
ylabel('LFP (mV)')
ax2 = nexttile;
h = pcolor(mt_t{n},mt_f,z_mt_s{n}); set(h,'EdgeColor','none'); colormap jet;
set(gca,'YScale','log','YTick',freqs(1:10:end),'CLim',[-1 5]), colorbar,
yline(bands.control{i, 1}(1), 'r');
yline(bands.control{i, 1}(2), 'r');
ylabel('Frequency (Hz)')
ax3 = nexttile;
hold on
plot(mt_t{n},bands.(field){i,2}); %bandpower over time
plot(mt_t{n},bands.(field){i,5}); %threshold
ylabel(['Bandpower ', bandlabel])
hold off
xlabel(tl, 'Time (s)');
title(tl, cond{n}, 'Interpreter', 'none');
linkaxes([ax1 ax2 ax3], 'x')
ax1.XLim = [0 lfp{2,n}(end)];
g = pan;
j = zoom;
set(g,'Motion','horizontal','Enable','on');
set(j,'Motion','horizontal','Enable','on');
%% reconstruct real times
for n = 1:size(lfp,2)
    driftperiods_all{1,n}(:,3) = driftperiods_all{1,n}(:,2)-driftperiods_all{1,n}(:,1);
end
prompt = 'What condition are you looking at? 1:control; 2:C21: 3:LFP_Ref: ';
n = input(prompt);
field = string(cond{n});
start = options.(field)(1); % control or C21
prompt = 'What time in the plot do you want to reconstruct in Labchart? time in seconds: ';
timeinplot = input(prompt);
a = start + timeinplot;
b = find(driftperiods_all{1,n}(:,1)<a, 1, 'last');
c = ceil(sum(driftperiods_all{1,n}(1:b,3)));
d = c + a;
while d > driftperiods_all{1,n}(b+1,1)
    b = find(driftperiods_all{1,n}(:,1)<d, 1, 'last');
    c = ceil(sum(driftperiods_all{1,n}(1:b,3)));
    d = c + a;
end
timeinlabchart = d
clear start timeinplot n a b c d prompt timeinlabchart
%% detection of oscillations using simple threshold
for n = 1:size(lfp,2) % control, C21 and LFP reference
    field = string(cond{n});
    for i = 1:size(bands.(field),1)
        bands.(field){i,6} = (bands.(field){i,2} > bands.(field){i,5});% thresholding bandpower
        bands.(field){i,7} = diff([0; bands.(field){i,6}]);% 1. derivative for onset and offset detection
        bands.(field){i,8} = find(bands.(field){i,7} == 1);% onsets (crossing threshold in positive direction)
        bands.(field){i,9} = find(bands.(field){i,7} == -1);% offsets (crossing threshold in negative direction)
        if length(bands.(field){i,8}) > length(bands.(field){i,9})
            bands.(field){i,8}(end) = [];
        end
        bands.(field){i,10} = bands.(field){i,9}-bands.(field){i,8};% length of events (how many data points above threshold)
        bands.(field){i,11} = bands.(field){i,10}(bands.(field){i,10} > MinEventTime*Spec_fs); % Minimum time an event is showing a power over threshold
        bands.(field){i,12} = bands.(field){i,8}(bands.(field){i,10} > MinEventTime*Spec_fs); % onsets of events, fullfiling the length criterium
        bands.(field){i,13} = bands.(field){i,9}(bands.(field){i,10} > MinEventTime*Spec_fs); % offsets of events, fullfiling the length criterium
        bands.(field){i,14} = floor(bands.(field){i,12}+bands.(field){i,11}./2); % event time defined as the middle of the event rounded towards zero
        bands.(field){i,15} = arrayfun(@(x)(find(mt_t{n}>(x./Spec_fs), 1, 'first')),bands.(field){i,14}); %how to get the correct time point?
        bands.(field){i,16} = mt_t{n}(bands.(field){i,15}+1)'; % event times in sec
        bands.(field){i,17} = diff(bands.(field){i,16}); % inter-event interval in sec, without drift
%         bands.(field){i,18} = goodidx{n}(bands.(field){i,15}); % event time index of original lfp before drift removal
%         lfp{3,n} = zeros(round(metadata.rectime{1, n}*metadata.lfp.fs), 1); % creating a vector of original lfp length
%         lfp{3,n}(bands.(field){i,18}) = 1; % translating event times back to original times
%         bands.(field){i,19} = find(lfp{3,n} == 1); % indices of real event times
%         bands.(field){i,20} = diff(bands.(field){i,19}); % inter-event interval as index, in real time, but longer drift removals might be problematic
    end
end

%% detection of oscillations using peakfinder
fr = 6; % fraction of 99th percentile(x0)-min(x0) for peak detection
for n = 1:size(lfp,2) % control, C21 and LFP reference
    field = string(cond{n});
    for i = 1:size(bands.(field),1)
        bands.(field){i,21} = smooth(bands.(field){i,2}, 8, 'sgolay',2); % smoothing bandpower
        % bands.(field){i,22} = peakfinder(bands.(field){i,21},(max(bands.(field){i,21})-min(bands.(field){i,21}))/fr,[],1,1,0); % defining sel with max-min divided by some factor
        bands.(field){i,22} = peakfinder(bands.(field){i,21}, diff(prctile(bands.(field){i,2}, [0 99]))/fr, std(bands.(field){i,2})*3,1,1,0); % defining sel with 99th percentile divided by some factor, with an additional threshold defined as 3xSD of bandpower
        
%         bands.(field){i,22} = peakfinder(bands.(field){i,21}, diff(prctile(bands.(field){i,2}, [0 99]))/fr,[],1,1,0); % defining sel with 99th percentile-min divided by some factor
        
        % bands.(field){i,22} = peakfinder(bands.(field){i,21}, std(bands.(field){i,2}),[],1,1,0); % defining sel with standard deviation
        bands.(field){i,23} = mt_t{n}(bands.(field){i,22}); % translating peak indices in times
        bands.(field){i,24} = diff(bands.(field){i,23}); % inter-event interval
    end
end

%% Extracting detected oscillatory events and calculate periodogram of each of them
win_osci = 1.5; % length of window to extract oscillations for event-based periodogram
for n = 1:(size(lfp,2)-1) % control, C21 without LFP_Ref
    field = string(cond{n});
    for i = 1:size(bands.(field),1)
        osci.(field){i,1} = int32(bands.(field){i,12} * win(2) * metadata.lfp.fs); % for onsets, extracting real detected events, begin of sequence
        osci.(field){i,1}(:,2) = int32(osci.(field){i,1}(:,1) + win_osci * metadata.lfp.fs - 1); % adding win_osci to onsets, end of sequence
        osci.(field){i,1}(:,3) = int32(bands.(field){i,13} * win(2) * metadata.lfp.fs); % for offsets, needed for using several criteria for event sorting
        if size(osci.(field){i,1},1) == 0 % if there is no event at least preallocate field in cell array
            osci.(field){i,2} = [];
            osci.(field){i,3} = [];
        else
            for m = 1:size(osci.(field){i,1},1)
                osci.(field){i,2}(:,m) = lfp{1,n}(osci.(field){i,1}(m,1):osci.(field){i,1}(m,2)); % extract lfp from detected events with a window length of win_osci
            end
        end
        [osci.(field){i,3}, win_osci_f] = periodogram(osci.(field){i,2}, [], [], metadata.lfp.fs); % calculate PSD of 1.5 s long LFP sequences of detected events only
        osci.(field){i,4} = mean(osci.(field){i,3},2); % calculate mean PSD over all detected events in that frequency range
        osci.(field){i,5} = std(osci.(field){i,3},0,2)/sqrt(size(osci.(field){i,3},2)); % calculate SEM of all PSDs from all detected events

        eventperiods = [];
        for s =1:size(osci.(field){i,1},1);
            eventperiods = [eventperiods, osci.(field){i,1}(s,1):osci.(field){i,1}(s,2)]; % vector of event indices, to be deleted from lfp
        end
        lfp{4,n} = lfp{1,n};
        lfp{4,n}(eventperiods) = []; % remove events from lfp
        osci.(field){i,6} = reshape(lfp{4,n}(1:end-rem(numel(lfp{4,n}),win_osci * metadata.lfp.fs)), win_osci * metadata.lfp.fs, []);  % generate table of 1.5 s lfp sequences outside of detected events
        [osci.(field){i,7}, win_osci_f] = periodogram(osci.(field){i,6}, [], [], metadata.lfp.fs); % calculate PSD of 1.5 s long LFP sequences of detected events only
        osci.(field){i,8} = mean(osci.(field){i,7},2); % same for sequences outside of events
        osci.(field){i,9} = std(osci.(field){i,7},0,2)/sqrt(size(osci.(field){i,7},2));        
    end
end
% dlmwrite([options.file_path, options.filename, '_', options.timestamp, convertStringsToChars(field), '.txt'], osci.(field){1,2}, '\t');
clear eventperiods m n i s 
%% SPW detection via LFP amplitude
for n = 1:size(lfp,2) % control, C21 and LFP reference
    field = string(cond{n});
    [SPW.(field){1,1}, SPW.(field){1,2}] = peakfinder(lfp{3,n}, (max(lfp{3,n})-min(lfp{3,n}))/30, SPW_thr, -1); % find LFP < SPW_thr
end
% sort out false positive detected SPWs
for n = 1:(size(lfp,2)-1) % control, C21 without LFP reference
    field = string(cond{n});
    eventperiods.(field) = [];
    % check if there is a relevant beta power increase
    for s =1:size(osci.(field){3,1},1); % only for beta band
        eventperiods.(field) = [eventperiods.(field), osci.(field){3,1}(s,1):osci.(field){3,1}(s,3)]; % vector of event indices from onset to offset, to be deleted from lfp
    end
    SPW.(field){3,1} = SPW.(field){1,1}(ismember(SPW.(field){1,1}, eventperiods.(field))); % keep only event indices, that are within event periods of elevated beta power
    SPW.(field){3,2} = SPW.(field){1,2}(ismember(SPW.(field){1,1}, eventperiods.(field))); % same for LFP amplitude
    % exclude close by SPWs with a distance smaller than 0.5 seconds
    a = SPW.(field){3,1}; % copy of SPW time
    a2 = SPW.(field){3,2}; % copy of SPW amplitude
    b = diff(a);
    while min(b) < metadata.lfp.fs/2 % as long as the minimum distance between SPWs is smaller than 0.5 seconds
        [c,d] = min(b); % c is the minimum distance and d is index of that value
        if a2(d)<=a2(d+1) % if amplitude of first SPW is larger than second
            a(d+1) = []; % the second event is deleted from the list of SPWs
            a2(d+1) = []; % same for LFP amplitude
        else 
            a(d) = [];
            a2(d) = [];
        end
        b = [diff(a)]; % regenerate diff vector of neighboring events
    end
    SPW.(field){4,1} = a;
    SPW.(field){4,2} = a2;
end
clear a a2 b c d s
%% this might be a try to detect SPWs, if the power increase in certain frequency bands doesn't work
% win_std = 0.075; % window length in sec for calculating the standard deviation of the raw lfp signal
% for n = 1:size(lfp,2)
%     field = string(cond{n});  
%     lfp_std{n} = movstd(lfp{1,n}, round(win_std*metadata.lfp.fs));%std over time of raw signal
% end
% LFP_Ref.lfp_std = movstd(LFP_Ref.lfp, round(win_std*metadata.lfp.fs));

%% Detecting HBOs
% for n = 1:size(lfp,2) % control, C21 and LFP reference
%     field = string(cond{n});
%     bands.(field){7,1} = '(4-8)-(10-16)';
%     bands.(field){7,2} = (bands.(field){1,2}-bp.LFP_Ref(1,1))-(bands.(field){6,2}-bp.LFP_Ref(6,1)); 
% end
% good try but doesn't work reliable 
%% Plot mean + SEM of all PSDs from events detected as power increases in various frequency bands
figure('Position',[300 100 1000 1000])
tl = tiledlayout('flow','TileSpacing','compact'); 
cond_legend = cond;
for i = 1:size(bands.(field),1)
    nexttile
    hold on
    for n = 1:(size(lfp,2)-1) % control, C21 without LFP reference
        field = string(cond{n});
        if isempty(osci.(field){i,4})
            continue
        else
        H = shadedErrorBar(win_osci_f, osci.(field){i,4}, osci.(field){i,5});
        set(H.mainLine, 'Color', cond_color(n,:))
        set(H.mainLine, 'LineWidth', 2)
        set(H.patch, 'FaceColor', cond_color(n,:))
        set(H.edge, 'Color', cond_color(n,:))
        cond_legend{n} = [cond{n}, ', n = ', num2str(size(osci.(field){i,3},2)), ' sequences'];
        end
    end
    legend(cond_legend, 'Interpreter', 'none')
    xlabel('Frequency [Hz]')
    ylabel('Power')
    xlim([3 40])
    bandlabel = ['[', num2str(bands.control{i, 1}(1)),'-',num2str(bands.control{i, 1}(2)), ' Hz]'];
    title(['Periodogram of events in ', bandlabel])
    hold off
end
title(tl, 'PSD of events detected by their power increase in various frequency ranges')
clear H 
export_fig(([options.file_path, options.filename, '_', options.timestamp]), '-pdf', '-silent', '-transparent', '-nocrop', '-append', '-noinvert') %-preserve_size
%% Plot mean + SEM of all PSDs from sequences that are not part of an event
figure('Position',[300 100 1000 1000])
tl = tiledlayout('flow','TileSpacing','compact'); 
for i = 1:size(bands.(field),1)
    nexttile
    hold on
    for n = 1:(size(lfp,2)-1) % control, C21 without LFP reference
        field = string(cond{n});
        if isempty(osci.(field){i,8})
            continue
        else
        H = shadedErrorBar(win_osci_f, osci.(field){i,8}, osci.(field){i,9});
        set(H.mainLine, 'Color', cond_color(n,:))
        set(H.mainLine, 'LineWidth', 2)
        set(H.patch, 'FaceColor', cond_color(n,:))
        set(H.edge, 'Color', cond_color(n,:))
        cond_legend{n} = [cond{n}, ', n = ', num2str(size(osci.(field){i,6},2)), ' sequences'];
        end
    end
    legend(cond_legend, 'Interpreter', 'none')
    xlabel('Frequency [Hz]')
    ylabel('Power')
    xlim([3 40])
    bandlabel = ['[', num2str(bands.control{i, 1}(1)),'-',num2str(bands.control{i, 1}(2)), ' Hz]'];
    title(['Periodogram of events in ', bandlabel])
    hold off
end
title(tl, 'PSD of event-free sequences')
clear H 
export_fig(([options.file_path, options.filename, '_', options.timestamp]), '-pdf', '-silent', '-transparent', '-nocrop', '-append', '-noinvert') %-preserve_size
%% Plotting LFP, spectrogram and detection results in bandpower
prompt = 'What condition do you want to plot? 1:control; 2:C21: 3:LFP_Ref: ';
n = input(prompt);
field = string(cond{n});
fprintf('1: theta %d-%d Hz\n',bands.control{1, 1})
fprintf('2: alpha %d-%d Hz\n',bands.control{2, 1})
fprintf('3: beta %d-%d Hz\n',bands.control{3, 1})
fprintf('4: theta-alpha %d-%d Hz\n',bands.control{4, 1})
fprintf('5: custom HBO %d-%d Hz\n',bands.control{5, 1})
fprintf('6: custom HBO %d-%d Hz\n',bands.control{6, 1})
prompt = 'What bandpower do you want to plot?: ';
i = input(prompt);
bandlabel = ['[', num2str(bands.control{i, 1}(1)),'-',num2str(bands.control{i, 1}(2)), ' Hz]'];
freqs = round(logspace(log10(1),log10(100),100),2);
gpu_lfp1 = gpuArray(lfp{1,n});
gpu_lfp2 = gpuArray(lfp{2,n});
gpu_lfp3 = gpuArray(lfp{3,n});
gpu_z_mt_s = gpuArray(z_mt_s{n});
q = figure('Position',[250 50 1200 900]); % [left bottom width height]
tl = tiledlayout(4,1,'TileSpacing','compact','Padding','none'); 
ax1 = nexttile;
hold on
plot(gpu_lfp2,gpu_lfp1);
plot(gpu_lfp2,medfilt1(lfp{1,n},metadata.lfp.fs));
plot(lfp{2,n}(SPW.(field){1,1}), SPW.(field){1,2}, 'o', 'MarkerEdgeColor', 'b')
hold off
ylabel('LFP (mV)')
ax2 = nexttile;
hold on
plot(gpu_lfp2,gpu_lfp3);
plot(lfp{2,n}(SPW.(field){4,1}), SPW.(field){4,2}, 'o', 'MarkerEdgeColor', 'r')
hold off
ylabel('LFP (mV)')
ax3 = nexttile;
h = pcolor(mt_t{n},mt_f,gpu_z_mt_s);
set(h,'EdgeColor','none'); colormap jet;
set(gca,'YScale','log','YTick',freqs(1:10:end),'CLim',[-1 5]), colorbar,
ylabel('Frequency (Hz)')
yline(bands.control{i, 1}(1), 'r');
yline(bands.control{i, 1}(2), 'r');
ax4 = nexttile;
hold on
plot(mt_t{n},bands.(field){i,2}); %bandpower over time
plot(mt_t{n},bands.(field){i,5}); %threshold
plot(bands.(field){i,16}, bands.(field){i,2}(bands.(field){i,14}), 'o', 'MarkerEdgeColor','r'); %detected events, with a larger bandpower than threshold
ylabel(['Bandpower ', bandlabel])
legend('raw bandpower', 'threshold', 'events using treshold')
hold off
ylim(ax4,'auto')
% ax5 = nexttile;
% hold on
% plot(mt_t{n},bands.(field){i,21}); %smoothed bandpower using Savitzky Golay method
% yline(std(bands.(field){i,2})*3);
% plot(bands.(field){i,23}, bands.(field){i,21}(bands.(field){i,22}), 'o', 'MarkerEdgeColor', 'r'); %detected events using peakfinder
% hold off
% ylabel(['Bandpower ', bandlabel])
% ylim(ax5,'auto')
% legend('smoothed bandpower', 'threshold', 'peakfinder')
% xlabel(tl, 'Time (s)');
title(tl, cond{n}, 'Interpreter', 'none');
linkaxes([ax1 ax2 ax3 ax4], 'x')
ax1.XLim = [0 lfp{2,n}(end)];
g = pan;
j = zoom;
set(g,'Motion','horizontal','Enable','on');
set(j,'Motion','horizontal','Enable','on');
savefig(q, [options.file_path, options.filename, '_', options.timestamp, '_', convertStringsToChars(field), '_', bandlabel], 'compact');

%% PSD of bandpower over time to detect oscillatory character
osci.win_welch = 60; % window length in sec for Welch's PSD estimate
osci.overlap = 0.5; % window overlap, number between 0 and 1
osci.win_length = round(osci.win_welch*Spec_fs); % larger windows increase resolution for smaller frequency ranges
osci.noverlap = round(osci.overlap*osci.win_length); % number of data points for window overlap
for n = 1:size(lfp,2) % control, C21 and LFP reference
    field = string(cond{n});
    for i = 1:size(bands.(field),1)
        [bands.(field){i,25}, osci.w{n}] = pwelch(bands.(field){i,2}, osci.win_length, osci.noverlap,[], Spec_fs, 'psd');
    end
end
%% Plotting PSD of bandpower over time
figure('Position',[300 100 1000 1000])
tl = tiledlayout('flow','TileSpacing','compact', 'padding', 'compact'); 
for i = 1:6 % size(bands.(field),1)
    nexttile
    hold on
    for n = 1:size(lfp,2)
        field = string(cond{n});
        plot(osci.w{n}, bands.(field){i,25}, 'color', cond_color(n,:))
    end
    legend(cond, 'Interpreter', 'none')
    xlabel('Frequency [Hz]')
    ylabel('Power')
    xlim([0.05 1])
    bandlabel = ['[', num2str(bands.control{i, 1}(1)),'-',num2str(bands.control{i, 1}(2)), ' Hz]'];
    title(['Bandpower ', bandlabel])
    hold off
end
title(tl, ['pwelch of bandpower over time (win = ', num2str(osci.win_welch), ' s) - oscillatory character'])
export_fig(([options.file_path, options.filename, '_', options.timestamp]), '-pdf', '-silent', '-transparent', '-nocrop', '-append', '-noinvert') %-preserve_size
%% Plotting a histogram of inter-event intervals, that is equivalent to the PSD of bandpower except that there are a lot more steps inbetween
figure('Position',[300 100 1000 1000])
tl = tiledlayout('flow','TileSpacing','compact', 'padding', 'compact'); 
for i = 1:6 % size(bands.(field),1)
    nexttile
    hold on
    for n = 1:(size(lfp,2)-1)
        field = string(cond{n});
        histogram(1./(bands.(field){i, 17}), 'Normalization', 'cdf', 'BinWidth', 0.005, 'DisplayStyle', 'stairs', 'EdgeColor', cond_color(n,:))
    end
    legend(cond, 'Interpreter', 'none')
    xlabel('Frequency [Hz]')
    ylabel('Count')
    xlim([0 1])
    bandlabel = ['[', num2str(bands.control{i, 1}(1)),'-',num2str(bands.control{i, 1}(2)), ' Hz]'];
    title(['Bandpower ', bandlabel])
    hold off
end
title(tl, ['Histogram of inter-event intervals - oscillatory character'])
export_fig(([options.file_path, options.filename, '_', options.timestamp]), '-pdf', '-silent', '-transparent', '-nocrop', '-append', '-noinvert') %-preserve_size
%% saving
tic
clear ax1 ax2 ax3 ax4 ax5 bandlabel field q g h i j m n serr prompt tl b cmap fr o p plot_color 
clear c data win_length data
clear e gpu_lfp1 gpu_lfp2 gpu_lfp3 gpu_z_mt_s
save([options.file_path, options.filename, '_', options.timestamp, '.mat'], '-v7.3', '-nocompression');
disp(['File was saved to ', ...
    options.file_path, options.filename, '_', options.timestamp, '.mat'])
toc
%% peek and see if the first few PCA or ICA bases give a hint that certain frequency bands explain variance across samples
[c s l] = pca(mt_s_log{2});
figure('Pos',[0 0 1200 500])
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact'); nexttile
plot(mt_f,s(:,1:3)'), title('PCA'), ylabel('raw mt spectra'), xlabel('Frequency (Hz)'), legend,
% nexttile
% rng default
% Mdl = rica(mt_s',6,'IterationLimit',1e4,'Standardize',false);
% plot(mt_f,Mdl.TransformWeights), title('ICA')
nexttile
[c s l] = pca(z_mt_s{2});
plot(mt_f,s(:,1:3)'),  ylabel('zscore mt spectra'), xlabel('Frequency (Hz)')


%% 
% Autocorrelation analysis to find discrete recurrene frequency of bursts
% 
% or power analysis on bandpower?

a = zeros(1, length(mt_t));
a(locs) = 1;
ac = dsp.Autocorrelator('Method', 'Time Domain');
y = ac(a);