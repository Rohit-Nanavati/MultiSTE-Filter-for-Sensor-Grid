function [newparticle,newweights,Dk,hopt] = resample_PF(particles,weights,D,regularizeMode)

    % % if nargin == 2
    % %     D = length(weights);
    % % end
    % 
    % CSW = cumsum(weights);
    % ii = 1 ;
    % indx = zeros(D,1);
    % u(1) = rand/D;
    % 
    % for jj = 1:D
    %     uj = u(1) + (jj-1)/D;
    %     while uj > CSW(ii)
    %         ii = ii + 1;
    %     end
    %     indx(jj) = ii;
    % end
    % 
    % fields = fieldnames(particles);
    % 
    % 
    % for findx = 1:length(fields)
    %     newparticle.(fields{findx}) = particles.(fields{findx})(indx,:);
    % end
    % % regularizeMode = 1;
    % 
    % % === ADD: carry mask through resampling ===
    % if isfield(particles,'mask')
    %     newparticle.mask = particles.mask(indx,:);
    %     % and immediately enforce Q=0 on inactive slots
    %     newparticle.Q = newparticle.Q .* (newparticle.mask ~= 0);
    % end


        % --- residual resampling ---
    % 1. expected counts
    expCount = D * weights(:);                  % D x 1

    % 2. deterministic copies
    copies_int = floor(expCount);               % integer part
    N_det = sum(copies_int);
    % indices for deterministic part
    det_idx = repelem((1:length(weights))', copies_int);  % N_det x 1

    % 3. stochastic copies for the remainder
    R = D - N_det;
    if R > 0
        % fractional remainder, renormalized
        frac = expCount - copies_int;
        if sum(frac) > 0
            frac = frac ./ sum(frac);
        else
            % all were integers somehow -> just uniform fallback
            frac = ones(size(frac)) / numel(frac);
        end

        % systematic draw on the remainder part
        CSW = cumsum(frac);
        ii = 1;
        stoch_idx = zeros(R,1);
        u0 = rand/R;
        for jj = 1:R
            uj = u0 + (jj-1)/R;
            while uj > CSW(ii)
                ii = ii + 1;
            end
            stoch_idx(jj) = ii;
        end

        full_idx = [det_idx; stoch_idx];        % total D x 1
    else
        full_idx = det_idx;
    end

    % now build newparticle from full_idx
    fields = fieldnames(particles);
    for f = 1:numel(fields)
        newparticle.(fields{f}) = particles.(fields{f})(full_idx,:);
    end

    % carry mask & enforce inactive Q=0 (keep your invariant)
    if isfield(newparticle,'mask')
        newparticle.Q = newparticle.Q .* (newparticle.mask ~= 0);
    end

    % equal weights after resample
    % newweights = ones(D,1) / D;



    [Np,N_sources] = size(particles.(fields{1}));

    if regularizeMode==1
        for ns = 1:N_sources
            xyparticlemat = [particles.x(:,ns) particles.y(:,ns)];
            Qparticlemat = particles.Q(:,ns);%*particles.mask(:,ns);
            xymeanest = weights'*xyparticlemat;
            Qmeanest = weights'*Qparticlemat;
            xyparticleCovariance = (xyparticlemat-xymeanest)'*diag(weights)*(xyparticlemat-xymeanest);
            QparticleCovariance = (Qparticlemat-Qmeanest)'*diag(weights)*(Qparticlemat-Qmeanest);
            Dk.xy(:,:,ns) = cholcov(xyparticleCovariance);
            Dk.Q(ns) = cholcov(QparticleCovariance);
        end
    elseif regularizeMode==2
        xyparticlemat = zeros(Np,2*N_sources);
        Qparticlemat = zeros(Np,N_sources);
        for ns = 1:N_sources
            xymatindx = 2*(ns-1)+1;
            xyparticlemat(:,xymatindx:xymatindx+1) = [particles.x(:,ns) particles.y(:,ns)];%.*particles.mask(:,ns);
            Qparticlemat(:,ns) = particles.Q(:,ns);%.*particles.mask(:,ns);
        end
        xymeanest = weights'*xyparticlemat;
        Qmeanest = weights'*Qparticlemat;
        xyparticleCovariance = (xyparticlemat-xymeanest)'*diag(weights)*(xyparticlemat-xymeanest);
        QparticleCovariance = (Qparticlemat-Qmeanest)'*diag(weights)*(Qparticlemat-Qmeanest);
        Dk.xy = cholcov(xyparticleCovariance);
        Dk.Q = cholcov(QparticleCovariance);
    else
        error("Regularization Mode Incorrectly set. Set 'regularizeMode' to 1 for source decouple regularisation and 2 for coupled regularisation");
    end

    newweights = ones(D,1)/D;

    mm = N_sources*length(fields); % The dimension of parameter space
    A = (4/(mm+2))^(1/(mm+4)); %InstantaneousGaussian
    hopt = A*(D^(-1/(mm+4)));

end