function [f,at1,at2,plotstruct] = mste_figplot_set(EstHistory,particles,tnow,f,at1,at2,plotstruct)

t_length = length(EstHistory);
measurementcell = cell(1,t_length);
N_sourcesmax = size(particles.Q,2);

for tindx = 1:t_length
    measurementcell{tindx} = EstHistory(tindx).measurement;
end

if ~isempty(EstHistory(t_length).scatterplotpoints)
    predN_sources = length(EstHistory(t_length).scatterplotpoints);
    partpoints = EstHistory(t_length).scatterplotpoints;
    estmat = EstHistory(t_length).confidentest.mean;
    for ns=1:N_sourcesmax
        if ns<=predN_sources
            set(plotstruct.histogramarray(ns),"Data",partpoints{ns}(:,3));
            set(plotstruct.histmean(ns),"Value",estmat(ns,3),'HandleVisibility','on');
            set(plotstruct.particle_plots(ns), 'XData', partpoints{ns}(:,1), 'YData', partpoints{ns}(:,2));
            set(plotstruct.sourestPlot(ns),'XData',estmat(ns,1),'YData',estmat(ns,2));
        else
            set(plotstruct.histogramarray(ns),"Data",nan);
            set(plotstruct.histmean(ns),"Value",nan,'HandleVisibility','off');
            set(plotstruct.particle_plots(ns), 'XData',nan, 'YData', nan);
            set(plotstruct.sourestPlot(ns),'XData',nan,'YData',nan);
        end
    end
end

timetxt = ['Time = ',num2str(tnow),'s and No. of Sources: ',num2str(length(EstHistory(end).estconfidence))];
title(at1,timetxt,'FontSize',12,'Interpreter','latex');
drawnow

    function [binedge,bincounts] = weightedHistogram(Wpnorm,theta,NoA)
        binedge = cell(NoA,1);
        bincounts = cell(NoA,1);
        for noaindx = 1:NoA
            maxval = max(theta(noaindx).Q);
            if maxval~=60
                maxval=60;
            end
            minval = min(theta(noaindx).Q);
            binwidth = 1;
            binedge{noaindx} = minval:binwidth:maxval;
            numBins = numel(binedge{noaindx}) - 1;
            bincounts{noaindx} = zeros(1,numBins);
            for ii = 1:numBins
                idx = theta(noaindx).Q >= binedge{noaindx}(ii) & theta(noaindx).Q < binedge{noaindx}(ii+1);
                bincounts{noaindx}(ii) = sum(Wpnorm(idx,noaindx));
            end
        end
    end
end