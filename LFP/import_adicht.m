function [f, data, comments] = import_adicht(file_path)
% Author: Jürgen Graf
% date: 23.05.2021
% dependent on:
% Installation Simple Data File SDK (provides a Windows dll)
%       found in C:\Program Files (x86)\ADInstruments\LabChart8\Extras
% ADInstruments (LabChart) SDK from Jim Hokanson
%       https://github.com/JimHokanson/adinstruments_sdk_matlab
%% Input
% file_path : str
%         Path of the file to read. An empty or missing input prompts the user.
%% Output
% f    
%         is an object that contains all metadata within the LabChart file
% data
%         cell array that contains raw data of selected channels during prompt input
%         data of different records are concatenated
% comments
%         table with comments noted during experiment
%         note that Time refers to the time in seconds since the start of the record
%% generate class adi.file

f = adi.readFile(file_path,'remove_empty_channels',false);

%% Generate table with channels to print in Command Window for channel selection
varNames = {'ChannelID', 'Name'};
varTypes = {'uint16', 'string'};
channels = table('Size', [length(f.channel_specs) 2], 'VariableTypes', varTypes, 'VariableNames', varNames);
clear varNames varTypes
for i = 1:length(f.channel_specs)
    channels(i,:) = {f.channel_specs(1, i).id, f.channel_specs(1, i).name};
end
%% Generate table with comments, noted during experiment
varNames = {'CommentID', 'recordID', 'Time', 'Comment'}; %Time in seconds since the start of the record
varTypes = {'uint16', 'uint16', 'single', 'string'};
n_comments = f.records(1, f.n_records).comments(1, end).id;
comments = table('Size', [n_comments 4], 'VariableTypes', varTypes, 'VariableNames', varNames);
clear varNames varTypes
n = 0; a = 0;
for r = 1:f.n_records
    if r == 1
        for i = 1:length(f.records(1,r).comments)
            n = n+1;
            comments(n,:) = {f.records(1,r).comments(1, i).id, f.records(1,r).comments(1, i).record, f.records(1,r).comments(1, i).time, f.records(1,r).comments(1, i).str};
        end
    elseif r>1
        a = seconds(datetime(f.records(1, r).data_start, 'ConvertFrom', 'datenum') - datetime(f.records(1, 1).data_start, 'ConvertFrom', 'datenum'));
        for i = 1:length(f.records(1,r).comments)
            n = n+1;
            comments(n,:) = {f.records(1,r).comments(1, i).id, f.records(1,r).comments(1, i).record, a + f.records(1,r).comments(1, i).time, f.records(1,r).comments(1, i).str};
        end
    end
end
clear n r i 
%% Select what channels to import
channels
prompt = 'What channels do you want to import? [1 3] for channel ID 1 and 3: ';
channelIDs = input(prompt);
clear prompt
%% That's the way suggested Jim Hokanson
% pres_chan = f.getChannelByName('Ephys');
% raw_pres_data = pres_chan.getData(1); % 1 stands for first RecordID
% fs = pres_chan.fs;
%% Import selected channels in a cell array
% channel objects are already existing in f, so I can directly use getData function
data = cell(1, f.n_channels);
for r = 1:f.n_records
    for i = channelIDs
        if isempty(data{1,i})
            data{1,i} = f.channel_specs(1, i).getData(r);
        else
            data{1,i} = [data{1,i}; f.channel_specs(1, i).getData(r)];
        end
    end
end
%% Message if n_records > 1
if f.n_records>1
    fprintf('Number of records is %i\n', f.n_records);
    fprintf('Records were concatenated at \n');
    for r = 1:f.n_records-1
        a = f.records(1, r).duration;
        b = seconds(datetime(f.records(1, r+1).data_start, 'ConvertFrom', 'datenum') - datetime(f.records(1, r).data_start, 'ConvertFrom', 'datenum')) - (f.records(1, r).duration);
        fprintf('%.2f seconds missing %.2f seconds\n', a, b);
    end
    fprintf('Keep in mind, that data are not continuous at this/these points \n');
end