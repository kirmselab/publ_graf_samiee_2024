%% This script runs Processing1_*.m for multiple directories.
%
%% Specify directories
% It is assumed that each of the following folders contains exactly N *_events.mat files, where
% N is the number experimental groups (e.g. 2 for a BEFORE vs. AFTER
% design).
folders ={
        '\\folder1\', ...
        '\\folder2\', ...
         };
%
outputDirectory = '\\outputfolder\'; % root output directory,
% where N new subfolders will be added (N, number of experimental groups)
outputSubdirectories = {
    'groupA_', ...
    'groupB_'}; % suffixes of subfolders corresponding to the N experimental groups
                  % (the sequence MUST reflect the sequence of *_events.mat files in 'folders')
%
%% Run Processing1 for each folder
timestamp = datestr(now, 'yymmdd_HHMMSS');
for k = 1:length(folders)
    clearvars -except folders k outputDirectory outputSubdirectories timestamp;
    Options.directory = folders{k};
    Processing1_Emx1GiDREADD_230906
end
%
%% Copy Processing1_*.mat files
disp('Copying Processing1_*.mat files to output directories ...')
for m = 1:length(outputSubdirectories)
    mkdir([outputDirectory, outputSubdirectories{m}, timestamp, '\']);
end
for k = 1:length(folders)
    filelist = dir([folders{k}, 'Processing1_*.mat']);
    if length(filelist) > length(outputSubdirectories) 
        error('The number of events.mat files does not match the number of experimental groups.')
    end
    CheckFinalROIs = cell(length(outputSubdirectories), 1);
    for m = 1:length(outputSubdirectories)
        aux = load([folders{k}, filelist(m).name], 'Events');
        % check if finalROIs in all 'Processing1_*.mat' files of a given
        % folder are equal
        CheckFinalROIs{m, 1} = aux.Events.finalROIs;
        if m > 1
            if ~isequal(CheckFinalROIs{m, 1}, CheckFinalROIs{1, 1})
                error('Final ROIs do not match.')
            end
        end
        copyfile([folders{k}, filelist(m).name], [outputDirectory, outputSubdirectories{m}, timestamp, '\']);
    end

end
disp('Completed.')
