function ibad = notSeegChannels(channel_labels,wm,varargin)
% function ibad = notSeegChannels(channel_labels,wm,bad_channels)
%   Returns indices of non-sEEG channels. It will look through
%   channel_labels and find all the channels that typically are not used in
%   our sEEG analyses (e.g. scalp EEG channels) and return indices for
%   these. If second input is set to 1, will also exclude channels with
%   white matter labels. If third input is given, it will also exclude
%   channels listed in bad_channels cell array.
%   
%   channel_labels (ch x [electrode name, anatomical label (optional)] cell array)
%   wm = 0 or 1 (include or exclude white matter channels?)
%   bad_channels = cell array of previously identified bad channel labels
%
%   DR 01/2023

ibad = [];
for ii = 1:length(channel_labels)
    elec = channel_labels{ii,1};
    if strncmp(elec,'C3',2) % exclude scalp EEG
        ibad = [ibad ii];
    elseif strncmp(elec,'C4',2)
        ibad = [ibad ii];
    elseif strncmp(elec,'Cz',2)
        ibad = [ibad ii];
    elseif strncmp(elec,'CZ',2)
        ibad = [ibad ii];
    elseif strncmp(elec,'Fp1',3)
        ibad = [ibad ii];
    elseif strncmp(elec,'Fp2',3)
        ibad = [ibad ii];
    elseif strncmp(elec,'F3',2)
        ibad = [ibad ii];
    elseif strncmp(elec,'F4',2)
        ibad = [ibad ii];
    elseif strncmp(elec,'F7',2)
        ibad = [ibad ii];
    elseif strncmp(elec,'F8',2)
        ibad = [ibad ii];
    elseif strncmp(elec,'Fz',2)
        ibad = [ibad ii];
    elseif strncmp(elec,'FZ',2)
        ibad = [ibad ii];
    elseif strncmp(elec,'P3',2)
        ibad = [ibad ii];
    elseif strncmp(elec,'P4',2)
        ibad = [ibad ii];
    elseif strncmp(elec,'Pz',2)
        ibad = [ibad ii];
    elseif strncmp(elec,'PZ',2)
        ibad = [ibad ii];
    elseif strncmp(elec,'O1',2)
        ibad = [ibad ii];
    elseif strncmp(elec,'O2',2)
        ibad = [ibad ii];
    elseif strncmp(elec,'T3',2)
        ibad = [ibad ii];
    elseif strncmp(elec,'T4',2)
        ibad = [ibad ii];
    elseif strncmp(elec,'T5',2)
        ibad = [ibad ii];
    elseif strncmp(elec,'T6',2)
        ibad = [ibad ii];
    elseif strncmp(elec,'LOC',3) % exclude EOG
        ibad = [ibad ii];
    elseif strncmp(elec,'ROC',3)
        ibad = [ibad ii];
    elseif strncmp(elec,'EMG',3) % exclude EMG
        ibad = [ibad ii];
    elseif strncmp(elec,'EKG',3) % exclude EKG
        ibad = [ibad ii];
    elseif strncmp(elec,'ECG',3)
        ibad = [ibad ii];
    elseif strncmp(elec,'SMELL',4) % exclude SMELL signal
        ibad = [ibad ii];
    elseif strncmp(elec,'RESP',4) % exclude RESP signal
        ibad = [ibad ii];
    elseif strncmp(elec,'DC',2) % exclude DC (analog input channels)
        ibad = [ibad ii];
    end
    if size(channel_labels,2)>1
        labl = channel_labels{ii,2};
        if isempty(labl) % exclude channels without label (usually outside of brain)
            ibad = [ibad ii];
        elseif contains(labl,'Unknown')
            ibad = [ibad ii];
        elseif wm && contains(labl,'White Matter') % exclude white matter-labeled channels
            ibad = [ibad ii];
        elseif wm && contains(labl,'White-Matter')
            ibad = [ibad ii];
        end
    end
end
if nargin==3
    bad_channels = varargin{1};
    bad = zeros(length(channel_labels),1);
    for ii = 1:length(bad_channels)
        bad = bad + contains(channel_labels(:,1),strtrim(bad_channels{ii}));
    end
    ibad = [ibad, find(bad)'];
end
ibad = sort(unique(ibad));