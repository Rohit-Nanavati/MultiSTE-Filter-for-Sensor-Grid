function particles = mergeCloseSourcesInParticles(particles, Domain, threshFrac, Qmin)
% MERGECLOSESOURCESINPARTICLES
% For each particle, find the closest pair among its K source locations.
% If distance < threshFrac * domainDiagonal, merge that pair:
%   - (x,y) := average; Q := sum
%   - re-seed the freed slot uniformly in Domain with Qmin
%
% Inputs:
%   particles : struct with fields x,y,Q of size [Np x K]
%   Domain    : polyshape defining the search area
%   threshFrac: scalar (e.g., 0.05 for 5% of domain diagonal)
%   Qmin      : scalar minimum release rate for re-seeded source
%
% Output:
%   particles : updated struct
    
indxwithzerosource = find(sum(particles.mask ~= 0,2)==0);
    Qth = Qmin;

    [Np, K] = size(particles.x);
    if K < 2, return; end

    % Domain diagonal (scale-free threshold)
    V    = Domain.Vertices;
    xmin = min(V(:,1)); xmax = max(V(:,1));
    ymin = min(V(:,2)); ymax = max(V(:,2));
    domDiag     = hypot(xmax - xmin, ymax - ymin);
    mergeThresh = threshFrac * domDiag;

    % Precompute mask for upper triangle (i<j) to exclude self distances
    U = triu(true(K), 1);             % KxK logical, diag=false, upper=true
    [Ui, Uj] = find(U);               % index lists for i<j pairs

    for i = 1:Np
        xi = particles.x(i, :).';
        yi = particles.y(i, :).';

        % Squared distances for all pairs (i<j) only
        dx = xi(Ui) - xi(Uj);
        dy = yi(Ui) - yi(Uj);
        d2 = dx.*dx + dy.*dy;

        % If K>=2, U is non-empty; handle safety anyway
        if isempty(d2), continue; end

        % Find minimal pair over i<j (self distances cannot appear here)
        [dmin2, idxMin] = min(d2);
        dmin = sqrt(dmin2);

        if dmin < mergeThresh
            aIdx = Ui(idxMin);
            bIdx = Uj(idxMin);

            % Merge positions (average) and release rates (sum)
            xnew = 0.5 * (particles.x(i,aIdx) + particles.x(i,bIdx));
            ynew = 0.5 * (particles.y(i,aIdx) + particles.y(i,bIdx));
            Qnew =        (particles.Q(i,aIdx) + particles.Q(i,bIdx));

            % Keep merged in aIdx; re-seed bIdx
            if logical(particles.mask(i,aIdx)) && logical(particles.mask(i,bIdx))
                particles.x(i,aIdx) = xnew; particles.y(i,aIdx) = ynew; particles.Q(i,aIdx) = Qnew;
    
                [xr, yr] = samplePointInPolyshape(Domain, xmin, xmax, ymin, ymax);
                particles.x(i,bIdx) = xr;   particles.y(i,bIdx) = yr;   particles.Q(i,bIdx) = Qmin;
                particles.mask(i,bIdx) = 0;
            else
                continue;
            end
        end
        if sum(particles.mask(i,:),2)>1 %particles with more than 1 active source
            oldmask = particles.mask(i,:);
            particles.mask(i,:) = double(oldmask & (particles.Q(i,:) > Qth));
            particles.Q(i,:)    = particles.Q(i,:) .* particles.mask(i,:);
        end
    end
    indxwithzerosource = find(sum(particles.mask ~= 0,2)==0);

    % % if Q is too small, mark its corresponding souce as inactive
    % Qth = 0.75;  % your threshold
    % if isfield(particles,'mask')
    %     indxgt1source = find(sum(particles.mask ~= 0,2)>1);
    %     modmask = (particles.mask ~= 0) & (particles.Q > Qth);
    %     particles.mask(indxgt1source,:) = modmask(indxgt1source,:);  % force inactive if Q<=Qth
    % else
    %     particles.mask = (particles.Q > Qth);                           % create mask if missing
    % end
    % 
    % indxwithzerosource = find(sum(particles.mask ~= 0,2)==0);
    % particles.mask = double(particles.mask);                            % keep 0/1 numeric
    % particles.Q    = particles.Q .* particles.mask;                     % ensure inactive sources contribute zero
    % 
    % indxwithzerosource = find(sum(particles.mask ~= 0,2)==0);

end

% -------- helper (same file) --------
function [xr, yr] = samplePointInPolyshape(Domain, xmin, xmax, ymin, ymax)
% Rejection sample uniformly in the polygon’s bounding box
    maxTries = 200;
    for t = 1:maxTries
        xr = xmin + (xmax - xmin) * rand;
        yr = ymin + (ymax - ymin) * rand;
        if isinterior(Domain, xr, yr)
            return;
        end
    end
    % Rare fallback: choose a random vertex
    V = Domain.Vertices;
    k = randi(size(V,1));
    xr = V(k,1); yr = V(k,2);
end
