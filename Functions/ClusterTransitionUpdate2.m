function [orderedparticles,clustermean] = ClusterTransitionUpdate2(particles)

[Np, N_sources] = size(particles.x);

Qmin = 0.75;                          % your threshold
% sel  = (particles.mask == 1) & (particles.Q > Qmin);   % Np×N_sources

xyQpoints = [particles.x(:), particles.y(:), particles.Q(:)*5];       % N_selected × 2

% xypoints = [particles.x(:) particles.y(:)];

[clusterID,clustermean] = kmeans(xyQpoints,N_sources);

clustermean = clustermean(:,end)/5;

permcombination = perms(1:N_sources);
permNumber = size(permcombination,1);

EstPosecostVector3d = zeros(permNumber,Np);
for m=1:permNumber
    partvec3d = zeros(N_sources,3,Np);
    partvec3d(:,1,:) = particles.x';
    partvec3d(:,2,:) = particles.y';
    partvec3d(:,3,:) = particles.Q';

    EstPosecostVector3d(m,:) = sum(vecnorm(partvec3d - clustermean(permcombination(m,:),:),2,2)); % maybe use mahalanobis distance?
end
[~,classificationIndx3d] = min(EstPosecostVector3d,[],1);

orderedparticles = particles;

if length(unique(classificationIndx3d))>1
    for np=1:Np
        orderedparticles.x(np,:) = particles.x(np,permcombination(classificationIndx3d(np),:));
        orderedparticles.y(np,:) = particles.y(np,permcombination(classificationIndx3d(np),:));
        orderedparticles.Q(np,:) = particles.Q(np,permcombination(classificationIndx3d(np),:));
        orderedparticles.mask(np,:) = particles.mask(np,permcombination(classificationIndx3d(np),:));
    end
end
end