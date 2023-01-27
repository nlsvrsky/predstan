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

% Suppressive Drive
idx = size(d,2);
if tau_t > 0
    temporalW = exp((-idx+1:0)*dt/tau_t) ./ sum(exp((-1e5:0)*dt/tau_t));
else
    temporalW = [zeros(1,idx-1) 1];
end
s = sum(sum(pool.*temporalW));

% Normalization
f = d(:,end) ./ (s + halfExp(sigma, p));

% Update firing rates
r = r_prev + (dt/tau)*(-r_prev + f);
