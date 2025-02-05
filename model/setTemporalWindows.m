function p = setTemporalWindows(p)

%% Sensory layer 1 (S1)
%  excitatory drive
if p.tauE1 > 0
    p.tempWE1 = exp(-(0:p.nt-1)*p.dt/p.tauE1) * (p.dt/p.tauE1);
else
    p.tempWE1 = [1 zeros(1,p.nt-1)];
end

%  suppresive drive
if p.tauS1 > 0
    p.tempWS1 = exp(-(0:p.nt-1)*p.dt/p.tauS1) * (p.dt/p.tauS1);
else
    p.tempWS1 = [1 zeros(1,p.nt-1)];
end

%% Sensory layer 2 (S2)
%  excitatory drive
if p.tauE2 > 0
    p.tempWE2 = exp(-(0:p.nt-1)*p.dt/p.tauE2) * (p.dt/p.tauE2);
else
    p.tempWE2 = [1 zeros(1,p.nt-1)];
end

%  suppresive drive
if p.tauS2 > 0
    p.tempWS2 = exp(-(0:p.nt-1)*p.dt/p.tauS2) * (p.dt/p.tauS2);
else
    p.tempWS2 = [1 zeros(1,p.nt-1)];
end

%% Voluntary attention layer
%  excitatory drive
if p.tauEAV > 0
    p.tempWEAV = exp(-(0:p.nt-1)*p.dt/p.tauEAV);
else
    p.tempWEAV = [1 zeros(1,p.nt-1)];
end

%  suppresive drive
if p.tauSAV > 0
    p.tempWSAV = exp(-(0:p.nt-1)*p.dt/p.tauSAV);
else
    p.tempWSAV = [1 zeros(1,p.nt-1)];
end

end