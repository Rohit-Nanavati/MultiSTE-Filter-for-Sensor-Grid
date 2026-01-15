function particles = mergeOrCullSources(particles, Domain, threshFrac)
% MERGEORCULLSOURCES
% For each particle (row):
%   - If the closest active pair is closer than threshFrac*domainDiagonal:
%       merge that pair (avg x,y; sum Q), then deactivate one of them.
%   - Else:
%       deactivate the active source with the smallest Q.
%
% Inputs:
%   particles : struct with fields x,y,Q,mask of size [Np x K]; mask==1 => active
%   Domain    : polyshape (used only to compute domain diagonal for threshold)
%   threshFrac: e.g., 0.05 for 5% of domain diagonal
%
% Output:
%   particles : updated (mask set to 0 when deactivated; Q set to 0)

[D, K] = size(particles.x);
if ~isfield(particles,'mask')
    error('particles.mask is required (1=active, 0=inactive).');
end
if K < 1, return; end

% Threshold uses domain diagonal
V         = Domain.Vertices;
domDiag   = hypot(max(V(:,1)) - min(V(:,1)), max(V(:,2)) - min(V(:,2)));
thresh2   = (threshFrac * domDiag)^2;   % compare squared distances (no sqrt)

for i = 1:D
    act = find(particles.mask(i,:) == 1);    % active column indices
    m = numel(act);
    if m < 1
        continue;  % nothing active
    elseif m == 1
        % Only one active source: by your rule “remove lowest Q” would zero out the last one.
        % Usually we keep at least one; skip to avoid going to zero.
        continue;
    else
        % Build pairwise squared distances among active columns (upper triangle)
        xi = particles.x(i, act).';
        yi = particles.y(i, act).';
        U = triu(true(m), 1);
        [u, v] = find(U);
        dx = xi(u) - xi(v);
        dy = yi(u) - yi(v);
        d2 = dx.*dx + dy.*dy;

        % Find closest pair
        [dmin2, idxMin] = min(d2);
        ai_local = u(idxMin); bi_local = v(idxMin);      % local (1..m)
        ai = act(ai_local);   bi = act(bi_local);        % global (1..K)
        if dmin2 < thresh2
            % ---- MERGE CLOSE PAIR ----
            xnew = 0.5 * (particles.x(i,ai) + particles.x(i,bi));
            ynew = 0.5 * (particles.y(i,ai) + particles.y(i,bi));
            Qnew =        (particles.Q(i,ai) + particles.Q(i,bi));
            particles.x(i,ai) = xnew;
            particles.y(i,ai) = ynew;
            particles.Q(i,ai) = max(Qnew, 0);           % keep non-negative
            particles.mask(i,bi) = 0;                   % deactivate the other
            particles.Q(i,bi)    = 0;                   % ensure no contribution
            % (positions of deactivated slot can be left as-is)
        else
            % ---- CULL LOWEST-Q ACTIVE ----
            [min_Q, localMinQ] = min(particles.Q(i, act));

            if min_Q < 2

                killCol = act(localMinQ);
                particles.mask(i, killCol) = 0;
                particles.Q(i, killCol)    = 0;
            else
                act_n = numel(act);
                randomIndex = randi(act_n);
                killCol = act(randomIndex);
                particles.mask(i, killCol) = 0;
                particles.Q(i, killCol)    = 0;

            end

        end
    end

end
end
