function [TransitParticles,TransitWeights,ProcessCov] = ParticleTransition(particles,weights,TransitionSigma,Domain,ts)

[Np, N_sources] = size(particles.x);
xmin = min(Domain.Vertices(:,1)); xmax = max(Domain.Vertices(:,1));
ymin = min(Domain.Vertices(:,2)); ymax = max(Domain.Vertices(:,2));

birthRate = 0.08;
deathRate = 0.08;

phi = ts.phi;
u = ts.u;


mask = (particles.mask ~= 0);
M = sum(mask, 2);                       % active count per particle
u = rand(Np,1);

doBirth = (u < birthRate) & (M < N_sources);         % only if an inactive slot exists
doDeath = (u > 1 - deathRate) & (M >= 2);  % else try death if ≥2 active

% ----- Birth: activate one inactive partition per selected particle -----
birthIdx = find(doBirth);
for ii = birthIdx.'
    inact = find(~mask(ii,:));          % inactive columns for this particle
    if ~isempty(inact)
        j = inact(randi(numel(inact))); % choose one to activate

        % Sample a new position inside the polygon (quick rejection in bbox)
        xr = xmin + (xmax - xmin)*0.05 + 0.9*(xmax - xmin)*rand; yr = ymin + (ymax - ymin)*0.05 +  0.9*(ymax - ymin)*rand;
        tries = 0;
        while ~isinterior(Domain, xr, yr) && tries < 10
            xr = xmin + (xmax - xmin)*0.05 + 0.9*(xmax - xmin)*rand; yr = ymin + (ymax - ymin)*0.05 +  0.9*(ymax - ymin)*rand;
            tries = tries + 1;
        end
        if ~isinterior(Domain, xr, yr)
            V = Domain.Vertices; k = randi(size(V,1)); xr = V(k,1); yr = V(k,2);
        end

        Q0 = gamrnd(2.5,6);

        particles.x(ii,j) = xr;
        particles.y(ii,j) = yr;
        particles.Q(ii,j) = Q0;
        mask(ii,j)        = true;
    end
end


% Write mask back
particles.mask = double(mask);          % keep numeric consistency if you rely on 0/1

% ----- Death: merge closest active pair and deactivate one slot -----
deathIdx = find(doDeath);

indxwithat2source = find(sum(particles.mask(deathIdx,:) ~= 0,2)<2);

if ~isempty(deathIdx)
    threshFrac = 0.15;  % e.g., 5% of domain diagonal

    % Build a small sub-struct for ONLY the particles that are “dying”
    sub.x    = particles.x(deathIdx, :);
    sub.y    = particles.y(deathIdx, :);
    sub.Q    = particles.Q(deathIdx, :);
    sub.mask = particles.mask(deathIdx, :);

    % Apply your rule: if a close pair exists -> merge; else -> cull lowest-Q
    sub = mergeOrCullSources(sub, Domain, threshFrac);

    % Write back
    particles.x(deathIdx, :) = sub.x;
    particles.y(deathIdx, :) = sub.y;
    particles.Q(deathIdx, :) = sub.Q;
    particles.mask(deathIdx, :) = sub.mask;
end
indxwithzerosource = find(sum(particles.mask ~= 0,2)==0);

%%%%
threshFrac = 0.1;   % 5% of domain diagonal
Qmin       = 0.75;   % minimal release rate for re-seeded source (tune as needed)

particles = mergeCloseSourcesInParticles(particles, Domain, threshFrac, Qmin);

%%%%%
% --------------------------
% Motion update (active only)
% --------------------------
TransitParticles = particles;           % start as copy

R_wg = [cos(phi) -sin(phi); sin(phi) cos(phi)];

for jj = 1:N_sources
    % Only update active rows
    activeRows = logical(particles.mask(:,jj));

    nx_w = TransitionSigma.x * randn(sum(activeRows),1);
    ny_w = TransitionSigma.y * randn(sum(activeRows),1);
    
    n_g = R_wg*[nx_w'; ny_w'];

    nx = n_g(1,:)';
    ny = n_g(2,:)';
    nq = TransitionSigma.Q * randn(sum(activeRows),1);

    % Propose
    TransitParticles.x(activeRows,jj) = particles.x(activeRows,jj) + nx;
    TransitParticles.y(activeRows,jj) = particles.y(activeRows,jj) + ny;
    TransitParticles.Q(activeRows,jj) = particles.Q(activeRows,jj) + nq;

    % Enforce domain & Q>0 for active rows (few retries)
    idx = activeRows & (TransitParticles.Q(:,jj) <= 0 | ...
                        TransitParticles.x(:,jj) < xmin | TransitParticles.x(:,jj) > xmax | ...
                        TransitParticles.y(:,jj) < ymin | TransitParticles.y(:,jj) > ymax);
    counter = 0;
    while any(idx) && counter < 5
        unfit = find(idx);

        nux_w = TransitionSigma.x * randn(numel(unfit),1);
        nuy_w = TransitionSigma.y * randn(numel(unfit),1);
        
        nu_g = R_wg*[nux_w'; nuy_w'];

        nux = nu_g(1,:)';
        nuy = nu_g(2,:)';

        nuq = TransitionSigma.Q * randn(numel(unfit),1);
        TransitParticles.x(unfit,jj) = particles.x(unfit,jj) + nux;
        TransitParticles.y(unfit,jj) = particles.y(unfit,jj) + nuy;
        TransitParticles.Q(unfit,jj) = particles.Q(unfit,jj) + nuq;

        idx = activeRows & (TransitParticles.Q(:,jj) <= 0 | ...
                            TransitParticles.x(:,jj) < xmin | TransitParticles.x(:,jj) > xmax | ...
                            TransitParticles.y(:,jj) < ymin | TransitParticles.y(:,jj) > ymax);
        counter = counter + 1;
    end

    % If still unfit, revert those rows to previous state (stable fallback)
    if any(idx)
        TransitParticles.x(idx,jj) = particles.x(idx,jj);
        TransitParticles.y(idx,jj) = particles.y(idx,jj);
        TransitParticles.Q(idx,jj) = particles.Q(idx,jj);
    end

    % Inactive rows: keep exactly as before (no motion)
    inactiveRows = ~activeRows;
    TransitParticles.x(inactiveRows,jj) = particles.x(inactiveRows,jj);
    TransitParticles.y(inactiveRows,jj) = particles.y(inactiveRows,jj);
    TransitParticles.Q(inactiveRows,jj) = particles.Q(inactiveRows,jj);
end

% --------------------------
% Likelihood process noise weight (unchanged)
% --------------------------
xparterror = TransitParticles.x - particles.x;
yparterror = TransitParticles.y - particles.y;
Qparterror = TransitParticles.Q - particles.Q;

particlesErrorMat = zeros(Np, 3*N_sources);
for ns = 1:N_sources
    xymatindx = 2*(ns-1)+1;
    particlesErrorMat(:,xymatindx)   = xparterror(:,ns);
    particlesErrorMat(:,xymatindx+1) = yparterror(:,ns);
end
particlesErrorMat(:,2*N_sources+1:3*N_sources) = Qparterror;

ProcessCov = blkdiag( kron(eye(N_sources), diag([TransitionSigma.x^2, TransitionSigma.y^2])), ...
                      TransitionSigma.Q^2 * eye(N_sources) );
TransitWeights = weights .* mvnpdf(particlesErrorMat, zeros(1,3*N_sources), ProcessCov);

end