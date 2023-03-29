function addSEEGlabels_rave
% function addSEEGlabels
%   Add brain segmentation labels to each sEEG contact and save in the 2nd
%   column of the channel_labels variable. Uses labels found in ITKSnap CSV
%   file.
%
%   DR 12/2022

% parameters
ddir = '/Users/garmo/Box/sEEG'; % data directory
subj = 'HUP253_RID864'; % subject

% load locations
cd(ddir);
locCSV = fullfile(subj,'electrodes.csv'); % CSV file path and name
try
    locations = readtable(locCSV);
catch
    error([locCSV ' not found']);
end

% match formatting with what is in channel_labels
locelec = locations.Label;
for ii = 1:length(locelec)
    blett = isstrprop(locelec{ii},'alpha');
    bnumb = isstrprop(locelec{ii},'digit');
    if sum(blett)>1 && sum(bnumb)==1 && ~strncmp(locelec{ii},'EMG',3) && ~strncmp(locelec{ii},'Fp',2)
        locelec{ii} = [locelec{ii}(blett) '0' locelec{ii}(bnumb)];
    end
end

% load channel labels
load(fullfile(subj,[subj '_Induction.mat']),'channel_labels'); % the CSV reconstructions are useful only for induction sessions (after EMU stay)

% add location to each channel label
for ii = 1:length(channel_labels)
    elec = channel_labels{ii,1};
    ind = find(strncmpi(locelec,elec,length(elec)));
    if ~isempty(ind)
        channel_labels{ii,2} = locations.FSLabel{ind};
    elseif any(ismember(elec,{'C3','C4','Cz','CZ','Fp1','Fp2','F3','F4','F7','F8','Fz','FZ','P3','P4','Pz','PZ','O1','O2','T3','T4','T5','T6'}))
        channel_labels{ii,2} = 'EEG';
    elseif any(ismember(elec,{'LOC','ROC'}))
        channel_labels{ii,2} = 'EOG';
    else
        channel_labels{ii,2} = '';
    end
end

% append to session .mat file
save(fullfile(subj,[subj '_Induction.mat']),'channel_labels','-append'); % overwrites if had run this function previously
