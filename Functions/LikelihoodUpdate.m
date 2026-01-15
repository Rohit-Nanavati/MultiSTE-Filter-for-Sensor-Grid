function [norm_weights,Likelihood] = LikelihoodUpdate(particles,weights,pos,measurement,D_threshold,ts,PDetect)

%% for a guassian likelihood function
Likelihood = ones(length(weights),1);

%% Minimal tweak: treat inactive sources as Q=0
if isfield(particles,'mask')
    psrc = particles;
    psrc.Q = particles.Q .* (particles.mask ~= 0);   % zero-out inactive sources
else
    psrc = particles;
end

%% original structure preserved
Likelihood = ones(length(weights),1);

for iter = 1:length(measurement)
    measurement_pred = Plume_model(psrc, pos(iter,:), ts);  % use masked-Q

    if measurement(iter) <= D_threshold
        NDsigma = D_threshold; % sigma for background noise
        likelihood = PDetect*(1/2)*(1 + erf((D_threshold - measurement_pred)./(NDsigma*sqrt(2)))) ...
                   + (1 - PDetect);
    else
        sigma = 0.3*measurement(iter) + D_threshold; % sigma for measurement when above threshold
        likelihood = (1./(sigma.*sqrt(2*pi))) .* exp(-((measurement_pred - measurement(iter)).^2) ./ (2*sigma.^2));
    end

    % light stability guard (no behavior change in normal cases)
    Likelihood = Likelihood .* max(likelihood, realmin);
end

weights = weights .* Likelihood + 10*eps;

% minimal normalization guard
s = sum(weights);
if s <= 0 || ~isfinite(s)
    norm_weights = ones(size(weights)) / numel(weights);
else
    norm_weights = weights ./ s;
end
end
