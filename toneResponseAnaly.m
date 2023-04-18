function toneResponseAnaly
% function toneResponseAnaly
%   Analysis of tone response data (run toneResponseSync.m first). Computes
%   time-continuous probability of response based on approach of Smith et
%   al. (2004)  Dynamic Analysis of Learning in Behavioral Experiments. The
%   Journal of Neuroscience 24(2):447-461.
%
%   DR 12/2022

% parameters
RTthresh = 2; % reaction time threshold for tone response (s)
        
% load tone response data
selpath = uigetdir([],'Select top directory of sEEG data');
cd(selpath);
datfile = dir(fullfile('**','*_Induction.mat'));
for ifile = 1:length(datfile)
    load(fullfile(datfile(ifile).folder,datfile(ifile).name),'annotations','t','tones');
    disp(datfile(ifile).name);
    if ~isstruct(tones)
        disp(tones)
        continue; % no tone response data
    end
    
    % compute reaction time (RT) to each tone (inf if no response)
    itrial = strcmp(annotations(:,1),'tone on');
    ttrial = cell2mat(annotations(itrial,2));
    if t(end)-ttrial(end) > 15 % fill in tone on events through end of session since task often stopped after LOC clearly obtained
        tint = [5 15]; % same inter-stimulus time interval range (s) used in toneResponse.m
        ct = ttrial(end);
        while ct < t(end)
            isi = tint(1) + diff(tint)*rand;
            ttrial(end+1) = ct + isi;
            ct = ct + isi;
        end
    end
    iresp = strcmp(annotations(:,1),'button press');
    tresp = cell2mat(annotations(iresp,2));
    RT = Inf*ones(1,length(ttrial)); % reaction time
    for jj = 1:length(ttrial)
        deltaR = tresp-ttrial(jj);
        deltaR(deltaR<=0) = [];
        deltaT = ttrial-ttrial(jj);
        deltaT(deltaT<=0) = [];
        if isempty(deltaR) || (deltaT(1)<deltaR(1)) % if no following button press or next tone occurs before next response, then no response
            continue;
        end
        RT(jj) = deltaR(1); % first button press after tone
    end
    
    % binarize responses based on RT threshold
    responses = RT<RTthresh;
    
    % binomial expectation maximization analysis (Smith et al 2004)
    runanalysis(responses,1,0.5);
    
    % loss of consciousness (LOC) criteria (Purdon et al 2013 PNAS)
    load('resultsindividual','pmid','p05','p95');
    pmid = pmid(2:end); p05 = p05(2:end); p95 = p95(2:end);
    ind = find(pmid < .05);
    if ~isempty(ind)
        LOC = ttrial(ind(1));
    else
        LOC = ttrial(end);
    end
    
    % save results in tones structure
    tones.RTthresh = RTthresh;
    tones.ttrial = ttrial;
    tones.responses = responses;
    tones.pmid = pmid; % response probability
    tones.p05 = p05; % response probability 95% confidence intervals
    tones.p95 = p95;
    tones.LOC = LOC;
    save(fullfile(datfile(ifile).folder,datfile(ifile).name),'tones','-append');

    % plot results
    figure('Name',mfilename,'NumberTitle','off','Units','normalized','Position',[1/4 1/4 1/2 1/2],'Color','w');
    patch([ttrial', fliplr(ttrial')],[p95 fliplr(p05)],'k','FaceColor',[0.7 0.7 0.7],'EdgeColor','none');
    hold on; plot(ttrial,pmid,'k','LineWidth',3);
    ind = find(responses);
    plot(ttrial(ind),ones(1,length(ind)),'bo','MarkerFaceColor','b','MarkerSize',8);
    ind = find(~responses);
    plot(ttrial(ind),zeros(1,length(ind)),'bo','MarkerFaceColor','b','MarkerSize',8);
    set(gca,'Box','off','TickDir','out','FontSize',14,'XLim',[ttrial(1) ttrial(end)],'YLim',[0 1]);
    line(get(gca,'XLim'),0.05*ones(1,2),'Color','k','LineStyle','--','LineWidth',2);
    line(LOC*ones(1,2),get(gca,'YLim'),'Color','r','LineStyle','-','LineWidth',2);
    text(LOC,max(get(gca,'YLim')),[' LOC = ' num2str(round(LOC)) ' s'],'HorizontalAlignment','left','VerticalAlignment','top','Color','r','FontSize',14);
    xlabel('time (s)')
    ylabel('response probability')
    title(datfile(ifile).name,'Interpreter','none');
    drawnow;
    orient landscape
    print toneResponseAnaly.ps -dpsc2 -fillpage -append
    if length(datfile) > 1
        pause(5); close;
    end
end
delete resultsindividual.mat
