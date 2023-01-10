function p = setStim(p)

% start and end times of T1
stimStart = p.stimOnset;
stimEnd = p.stimOnset + p.stimDur;

% make stim
timeSeries = zeros([p.norient p.nt]);
timeSeries(p.stimseq(1),unique(round((stimStart:p.dt:stimEnd)/p.dt))) = 1; % T1
timeSeries(p.stimseq(2),unique(round(((stimStart:p.dt:stimEnd) + p.soa)/p.dt))) = 1; % T2

% stimulus contrasts
contrSeries = zeros([p.norient p.nt]);
contrSeries(p.stimseq(1),unique(round((stimStart:p.dt:stimEnd)/p.dt))) = p.contrast(1); % T1
contrSeries(p.stimseq(2),unique(round(((stimStart:p.dt:stimEnd) + p.soa)/p.dt))) = p.contrast(2); % T2

p.stim = timeSeries;
p.stimContrast = contrSeries;