function p = setStim(p)

if strcmp(p.stimMode,'standard')
    % start and end times of T1
    stimStart = p.stimOnset;
    stimEnd = p.stimOnset + p.stimDur; 

    % cue present?
    timeSeries = zeros([p.norient p.nt]); % use the orientations as cues for now
    if p.stimseq(1) > 0 
        timeSeries(p.stimseq(1),unique(round((stimStart:p.dt:stimEnd)/p.dt))) = 1; % T1
    end

    % reward present?
    rwdSeries = zeros([1 p.nt]); % only one reward type for now
    if p.stimseq(2) > 0
        rwdSeries(p.stimseq(2),unique(round(((stimStart:p.dt:stimEnd) + p.soa)/p.dt))) = 1;
    end

    % stimulus contrasts
    contrSeries = zeros([p.norient p.nt]);
    if p.stimseq(1) > 0
        contrSeries(p.stimseq(1),unique(round((stimStart:p.dt:stimEnd)/p.dt))) = p.contrast(1); % T1
    end
    %contrSeries(p.stimseq(2),unique(round(((stimStart:p.dt:stimEnd) + p.soa)/p.dt))) = p.contrast(2); % T2

    p.stim = timeSeries;
    p.stimContrast = contrSeries;
    p.rwd = rwdSeries;
else % random stim array
    p.stim = zeros(4,length(p.tlist));
    p.stim(1,:) = randi([0 1],1,length(p.tlist));
    p.stimContrast = p.stim;
end
