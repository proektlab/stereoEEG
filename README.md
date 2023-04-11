# stereoEEG
Repository to house code for a collaborative project analyzing stereoEEG data during anesthesia induction, anesthesia emergence, and sleep for NINDS-R01NS113366: [The role of dynamical criticality in human perception](https://reporter.nih.gov/project-details/9969278)

## Downloading data from iEEG.org
Process and save .mat file for each session with:  
`process_ieeg_data;`  
If tone response data is available for the session, make sure the 'toneResponse_Induction' subfolder is in the specified 'data_dir/subject' path so that it is included in the .mat file.

## Add brain segmentation labels to each sEEG contact
Run `addSEEGlabels;` (or `addSEEGlabels_rave;` for newer subjects) to add labels to the channel_labels variable in the .mat file. Uses labels found in ITKSnap CSV file (addSEEGlabels) or exported from RAVE (addSEEGlabels_rave). This is currently only relevant for the induction sessions.

## Sync behavioral-sEEG data
If tone response data is available for the session, run `toneResponseSync;`. 
This adds the tone on, tone off, button press, button release timestamped events to the session annotations. 

## Analyze behavioral data
Run `toneResponseAnaly;` (requires functions in Smith2004 subfolder) to estimate probability of response based on response data following analysis of Smith et al 2004 J Neurosci. Estimate of time of loss of consciousness (LOC) follows similar criteria as Purdon et al 2013 PNAS.

## Identify channels to include in analyses
Run `ibad = notSEEGchannels(channel_labels,wm,bad_channels);` to return indices of channels not to include in subsequent analyses. Removes non-sEEG channels (e.g. scalp EEG channels), channels labeled as white matter if set second argument to 1, and any previously identified 'bad' channels listed in cell array in third argument.

## Qualitative analysis: time domain
Load .mat file for session. Then produce scrolling time plot with:  
`timePlot(data,sample_rate,[],[],1,annotations,channel_labels(:,1));`. 

## Qualitative analysis: frequency domain
Load .mat file for session. Then produce power spectral density (PSD) plot with:  
`pwelchPlot(session.data,session.sample_rate,[0.5 500],5,1,500,session.channel_labels(:,1));`  
Identify bad channels by clicking on outlier PSD with Edit Plot tool and finding channel name:  
`get(gco,'Tag')`

## Stability analysis 
Replicating analysis of Solovey et al 2015 J Neurosci (Figures 2, 3, 4).
Load .mat file for session. Then perform analysis with:  
`seegStability(session);`.
All available sessions can be batch processed using:  
`seegStabilityBatch;` which will generate a plot for each session and save a .mat file for population analysis.  
After loading that population .mat file containing the structure `stability`, the population analysis can be run with:  
`seegStabilityPop(stability);`

## Time-frequency analysis
Computes both average power and average criticality in each time x frequency bin. 
Load .mat file for session. Then perform analysis with:  
`seegTF(session);`.

## Criticality-frequency analysis
Computes distribution of mode criticality vs mode frequency for a specified pair of time windows (awake vs anesthetized) in a session. 
Load .mat file for session. Then peform analyis with:
`seegCF(session,tawake,tanesth);`.
