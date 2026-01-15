function [confidentest,estconfidence,totalest,totalestconfidence,clustermean,confidentrelparticles] = EstimateCompute(particles,weights)
    
    [particles,clustermean] = ClusterTransitionUpdate2(particles);
    
    n_sourcemax = size(particles.x,2);
    totalest = zeros(n_sourcemax,3);
    totalestconfidence = zeros(n_sourcemax,1);
    totalestuncertainity = zeros(n_sourcemax,3);
    relparticles = cell(n_sourcemax,1);

    for ns = 1:n_sourcemax
        activeindx = particles.mask(:,ns)~=0;
        relweights = weights(activeindx);
        relparticlesmat = [particles.x(activeindx,ns) particles.y(activeindx,ns) particles.Q(activeindx,ns)];
        normrelweights = relweights/sum(relweights);
        totalest(ns,:)  = normrelweights'*relparticlesmat;
        totalestuncertainity(ns,:) = normrelweights'*((relparticlesmat - totalest(ns,:)).^2);
        totalestconfidence(ns) = sum(relweights)*100;
        
        relparticles{ns} = relparticlesmat;
    end
    
    confidentEstindx = totalestconfidence>50;
    confidentest.mean = totalest(confidentEstindx,:);
    confidentest.var = totalestuncertainity(confidentEstindx,:);
    confidentest.normedstd = sqrt(sum(confidentest.var,'all')/sum(confidentEstindx));
    confidentest.normedPosestd = sqrt(sum(confidentest.var(:,1:2),'all')/sum(confidentEstindx));
    estconfidence = totalestconfidence(confidentEstindx);
    confidentrelparticles = relparticles(confidentEstindx);
end