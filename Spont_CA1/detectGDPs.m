function [GDPIndex, GDP_Onsets, GDP_Offsets, data_sorted, options] = detectGDPs(data, Nrois, options)
% This functions detects GDPs from Ca2+ transient (CaT) times as obtained from multi-neuronal Ca2+ imaging.
% ==========================================================================================================
% Authors:
% Chuanqiang Zhang, Knut Kirmse (2022), modified 09/2023
% ==========================================================================================================
% Reference:
% The definition follows Flossmann et al. Somatostatin Interneurons Promote
% Neuronal Synchrony in the Neonatal Hippocampus. Cell Rep (2019)
% --
% "Network events were operationally deﬁned as GDPs as follows: (1) CaTs were classiﬁed as GDP-related if 
% they fell into a 500-ms-long time interval during which the fraction of active cells (i.e., cells with
% >= 1 detected CaT) was >= 20%. (2) Neighboring GDP-related CaTs were assigned to the same GDP if both
% shared a 500-ms-long time interval during which the fraction of active cells was >= 20%,
% otherwise to separate GDPs."
% ==========================================================================================================
% Input:
% • data: a nx2 array with
%   ... CaT times [ms] in col1 and
%   ... corresponding ROI numbers [] in col2
% • Nrois: a scalar specifying the total number of analyzed ROIs
% • options: a struct with the following mandatory fields
%   ... options.interval                     % time window [ms]
%   ... options.threshold_activecells        % threshold fraction of active cells []
% ==========================================================================================================
% Output:
% • GDPIndex: a vector (same length as data(:, 1)), in which GDPIndex(m) specifies the running
%   number of the GDP to which the m-th CaT (in data) belongs. If the
%   m-th CaT does not belong to any GDP, then GDPIndex(m) = NaN.
% • GDP_Onsets: vector of GDP onsets [ms]
% • GDP_Offsets: vector of GDP offsets [ms]
% • data_sorted: data (input) but sorted by col1 then col2
% • options: see above, with the following additonal field
%   ... options.threshold_count              % threshold number of active cells []
% ==========================================================================================================
%
%% Step 1: Determine if each CaT belongs to any GDP
options.threshold_count = options.threshold_activecells * Nrois;
data_sorted = sortrows(data, [1, 2]);
Nevents = size(data_sorted, 1);
GDPflag_matrix = NaN(Nevents, Nevents);
for k = 1:Nevents
    search_onset = data_sorted(k, 1);
    search_offset = data_sorted(k, 1) + options.interval;
    EventInInterval = data_sorted(:, 1) >= search_onset & data_sorted(:, 1) < search_offset;  
    aux = data_sorted(:, 2);
    cellcount = numel(unique(aux(EventInInterval)));
    if cellcount >= options.threshold_count
        GDPflag_matrix(:, k) = EventInInterval;
    else
        GDPflag_matrix(:, k) = false(Nevents, 1);
    end
end
GDPflag = max(GDPflag_matrix, [], 2);
%
%% Step 2: Determine onsets and offsets of individual GDPs (KK method)
a = movmin(GDPflag_matrix, [0, 1], 1);
a = max(a, [], 2);
a = [diff([0; a])];
GDP_Onsets = find(a == 1);
GDP_Offsets = find(a == -1);
NBursts = numel(GDP_Onsets);
GDPIndex = NaN(Nevents, 1); % pre-allocate
if NBursts > 0
    if numel(GDP_Offsets) < NBursts
        % if the last CaT belongs to a GDP, this event is defined as a GDP offset
        GDP_Offsets = [GDP_Offsets; Nevents];
    end
    for k = 1:NBursts
        GDPIndex(GDP_Onsets(k):GDP_Offsets(k), 1) = k;
    end
    GDP_Onsets = data_sorted(GDP_Onsets, 1); % [ms]
    GDP_Offsets = data_sorted(GDP_Offsets, 1); % [ms]
end
%
end