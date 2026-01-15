function [updated_particles,updated_weights] = PFupdateSimple(particlesMinus,weightsMinus,measurement,pos,Np,D_threshold,ts,PDetect,TransitionSigma,Domain,clusterstatus)


% Predicted Measurements at the current pos as per the predicted source
... location from the particle distribution

regularizeMode = 2;

%==========================================================================

% --- Merge close sources per particle (BEFORE process model) ---
% threshFrac = 0.05;   % 5% of domain diagonal
% Qmin       = 0.01;   % minimal release rate for re-seeded source (tune as needed)
% 
% particlesMinus = mergeCloseSourcesInParticles(particlesMinus, Domain, threshFrac, Qmin);
%==========================================================================

[particles,~,~] = ParticleTransition(particlesMinus,weightsMinus,TransitionSigma,Domain,ts);


[weights,~] = LikelihoodUpdate(particles,weightsMinus,pos,measurement,D_threshold,ts,PDetect);


Neff = 1/sum(weights.^2);

if Neff < 0.6*Np
%     disp("Particles Resampled");
    [particles,weights,~,~] = resample_PF(particles,weights,Np,regularizeMode);

    % === ADD: enforce Q=0 after resample ===
    if isfield(particles,'mask')
        particles.Q = particles.Q .* (particles.mask ~= 0);
    end
    
end 

updated_particles = particles;
updated_weights = weights;


end