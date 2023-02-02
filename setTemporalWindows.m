function p = setTemporalWindows(p)

%% Sensory layer 1 (S1)
%  excitatory drive
if p.tauE1 > 0
    tempW = exp(-(1:p.nt)*p.dt/p.tauE1);
    p.tempWE1 = tempW ./ sum(tempW);
else
    p.tempWE1 = [1 zeros(1,p.nt-1)];
end

%  suppresive drive
if p.tauS1 > 0
    tempW = exp(-(1:p.nt)*p.dt/p.tauS1);
    p.tempWS1 = tempW ./ sum(tempW);
else
    p.tempWS1 = [1 zeros(1,p.nt-1)];
end

%% Sensory layer 2 (S2)
%  excitatory drive
if p.tauE2 > 0
    tempW = exp(-(1:p.nt)*p.dt/p.tauE2);
    p.tempWE2 = tempW ./ sum(tempW);
else
    p.tempWE2 = [1 zeros(1,p.nt-1)];
end

%  suppresive drive
if p.tauS2 > 0
    tempW = exp(-(1:p.nt)*p.dt/p.tauS2);
    p.tempWS2 = tempW ./ sum(tempW);
else
    p.tempWS2 = [1 zeros(1,p.nt-1)];
end

end