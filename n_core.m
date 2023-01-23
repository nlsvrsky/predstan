function [r, f, s] = n_core(d, sigma, p, r_prev, tau, tau_t, dt)

% function [r, f, s] = n_core(d, sigma, p, r_prev, tau, dt)
%
% r_prev = p.r(:,idx-1)
% d = d
% s = p.s(:,idx)
% f = p.f(:,idx)
% r = p.r(:,idx)

% We let the suppressive drive be broad in feature space by letting all
% features contribute equally to the normalization pool

% Normalization pool
pool = abs(d); % abs in case drive is negative

% Excitatory drive
idx = size(d,2)-1;
if tau_t > 0
    e = sum(pool .* exp((-idx:0)*dt/tau_t),2) ./ sum(exp((-1e5:0)*dt/tau_t));
else
    e = pool(:,end);
end

% Suppressive Drive
s = sum(abs(e));

% Normalization
f = e ./ (s + halfExp(sigma, p));

% Update firing rates
r = r_prev + (dt/tau)*(-r_prev + f);
