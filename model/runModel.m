function [perfv, p, ev] = runModel(opt, modelClass, rsoa, rseq, rcond)

% function [perfv, p, ev] = runModel(opt, modelClass, rsoa, rseq, rcond)
%
% INPUTS
% opt
%   A structure with parameter values. Any specified values will
%   overwrite the defaults in setParameters. Also use opt to pass in
%   display options for runModel.
% modelClass 
%   The model name, see setParameters
% rsoa 
%   SOA conditions to run, options listed in runModel
% rseq
%   Orientation sequence conditions to run, options listed in runModel
% rcond
%   Attentional cueing conditions to run, options listed in runModel
%
% All input arguments are optional.
%
% OUTPUTS
% perfv
%   Performance (decision evidence) for each target (in cells), cueing
%   condition, and SOA. Organized for plotting.
% p
%   Structure containing parameter values, see setParameters
% ev
%   Decision evidence for each condition combination

%% Set opt and modelClass
if nargin < 1
    opt = [];
end
if nargin < 2 || isempty(modelClass)
    modelClass = '';
end

%% Set params
p = setParameters(opt, modelClass);

%% Display
if isfield(opt,'display')
    display = opt.display;
    opt = rmfield(opt,'display');
    
    if isfield(display,'plotTS')
        plotTS      = display.plotTS;
    end
    if isfield(display,'plotPerf')
        plotPerf    = display.plotPerf;
    end
    if isfield(display,'verbose')
        verbose = display.verbose;
    end
end

if ~exist('plotTS','var') 
    plotTS          = 0;
end
if ~exist('plotPerf','var')
    plotPerf        = 1;
end
if ~exist('verbose','var')
    verbose         = 0;
end

%% Set conditions to simulate
condnames       =  {'cueT1','cueT2','cueN'};

% Stimuli
tilt = 2*pi/180;
orientations = [pi/2+tilt pi/2-tilt tilt pi-tilt]; % V-CCW V-CW H-CCW H-CW
p.norient = numel(orientations);

% Conditions
contrasts = p.stimContrasts;
%soas      = [100:50:500 800]; % why was this commented out?
stimseqs  = {[2 0], [2 1], [0 1]}; %[cue reward] %{[2 1],[2 2],[2 3],[2 4]};

% Pick conditions to run
rcontrast = 1:size(contrasts,2); % contrast levels to run
ncontrast = numel(rcontrast);

if ~exist('rcond','var') || isempty(rcond)
    rcond = 1:3; % conditions to run
end
ncond = numel(rcond);

if ~exist('rsoa','var') || isempty(rsoa)
    %rsoa = 1:numel(soas); % soa levels to run COMMENTED OUT WHY?
    rsoa = 250;
end
nsoa = numel(rsoa);

if ~exist('rseq','var') || isempty(rseq)
    rseq = 4; % sequences to run
end
nseq = numel(rseq);

% RF response
for iO = 1:p.norient
    p.rfresp(iO,:) = rfResponse(orientations(iO), p.ntheta);
end

%% loop through all conditions to run
ev = zeros(2,nsoa,ncond,ncontrast);
for icond = 1:numel(rcond)
    cond = rcond(icond);
    condname = condnames{cond};
    
    %% Loop through contrast levels
    for icontrast = 1:ncontrast
        c = rcontrast(icontrast);
        for isoa = 1:nsoa
            s = rsoa(isoa);
            for iseq = 1:nseq
                % set conditions
                q = rseq(iseq);
                p.cond = cond;
                p.condname = condname;
                p.contrast = contrasts(:,c);
                p.soa = s; %soas(:,s);
                p.nstim = numel(p.soa)+1;
                p.stimseq = stimseqs{q};
                if p.stimseq > 0
                    p.orientseq = orientations(p.stimseq);
                end
                if verbose
                    fprintf('cond: %s contrast: %1.2f soa: %d seq: %d %d\n\n', condname, p.contrast, p.soa, p.stimseq)
                end
                
                % distribute voluntary attention
                if p.distributeVoluntary
                    if strcmp(condname, 'cueN')
                        w = p.AVNeutralT1Weight;
                    else
                        w = p.AVProp;
                    end
                    p.AVWeights = distributeAttention(p.soa, p.tR, w);
                end
                
                % set time series
                p = initTimeSeries(p);
                p = setStim(p);
                p = setTask(p,condname);
                p = setDecisionWindows(p);
                p = setTemporalWindows(p);

                % run the model
                p = n_model(p);

                % get evidence
                for iStim = 1:2
                    p.ev(:,iStim) = p.rd(iStim,end);
                end
                
                % make evidence positive for the correct stimulus
                p.ev = p.ev.*(-1).^p.stimseq; % flip sign for CCW (odd #s)
                
                % scale evidence (for fitting)
                p.ev = p.ev.*[p.scaling1 p.scaling2];
                
                % store evidence
                ev(:,isoa,icond,icontrast,iseq) = p.ev;
                
                %% Draw time series
                if plotTS
                    plotTimeSeries(p, condname)
                end
            end
        end
    end
end

%% plot multiple conditions
if plotPerf
    for icontrast = 1:numel(rcontrast)
        perfv = plotPerformance(condnames(rcond), rsoa, mean(ev(:,:,:,icontrast,:),5)); % soas(rsoa) -> rsoa
    end
else
    perfv = [];
end
