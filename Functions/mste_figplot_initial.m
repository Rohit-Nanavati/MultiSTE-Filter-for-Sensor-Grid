function [f,at1,at2,plotstruct] = mste_figplot_initial(EstHistory,Domain,ts,fighandle,tnow,N_sources)

xmin = min(Domain.Vertices(:,1));
xmax = max(Domain.Vertices(:,1));
ymin = min(Domain.Vertices(:,2));
ymax = max(Domain.Vertices(:,2));
t_length = length(EstHistory);
NoA = size(EstHistory(1).pose,1);

measurementcell = cell(1,t_length);
hmap = brewermap(2*NoA+1,'Set1');
hmap = hmap([1:5,7:end],:);

f = figure(fighandle);
f.Units = 'normalized';
f.OuterPosition = [0.302604166666667,0.1525,0.334375,0.820833333333333]; % 0.348958333333333,0.034259259259259,0.646875,0.82037037037037
pose_vec = EstHistory(1).pose;
for tindx = 1:t_length
    measurementcell{tindx} = EstHistory(tindx).measurement;
end

at1 =  subplot(3,1,1); hold on;% grid on;
at1.Position = [0.148402555910543,0.771684981684982,0.75,0.16]; %0.148402555910543,0.771684981684982,0.75,0.16
plotstruct.histogramarray = [];
plotstruct.histmean = [];
for ns=1:N_sources
    h = histogram(at1,nan,'Normalization','pdf','DisplayStyle','bar','HandleVisibility','off','EdgeColor','none','FaceColor',hmap(NoA+ns,:),'FaceAlpha',0.3,'BinLimits', [0, 20]);
    plotstruct.histogramarray = [plotstruct.histogramarray;h];
    plotstruct.histmean = [plotstruct.histmean; xline(NaN,'--','Color',hmap(NoA+ns,:),'linewidth',2,'Alpha',1)]; %,'LabelHorizontalAlignment','center','Interpreter','latex','LabelVerticalAlignment','bottom');
    legendlabelsQ{ns} = ['$\hat{\Theta}(Q)_{',num2str(ns),'}$'];
end
xline(ts.Q,'-','Color','k','linewidth',2,'Alpha',0.5);
xlabel('$Q$'); ylabel('$p(Q)$');
legendlabelsQ{end+1} = 'Truth';
legend(legendlabelsQ,'Interpreter','latex','Orientation','horizontal','FontSize',12,'Location','bestoutside',Position=[0.127385796412274,0.65914797008547,0.768522350993377,0.03710699023199]);

at2 =  subplot(3,1,[2 3]); hold on;
at2.Position = [0.151090342679128,0.09406779661017,0.72702492211838,0.527401129943503];

plot(pose_vec(:,1),pose_vec(:,2),'Marker','x','Color','k','LineStyle','none');
legendlabelsS = {'Sensor Network'};

plotstruct.particle_plots=[];
for ns=1:N_sources
    plotstruct.particle_plots = [plotstruct.particle_plots;scatter(NaN,NaN,0.5,'o','HandleVisibility','off','MarkerEdgeAlpha',0.15,'MarkerEdgeColor',hmap(NoA+ns,:),'MarkerFaceAlpha',0.25,'MarkerFaceColor',hmap(NoA+ns,:))]; hold on;
end

plotstruct.sourestPlot = [];
for ns=1:N_sources
    plotstruct.sourestPlot = [plotstruct.sourestPlot; plot(nan,nan,'p','MarkerSize',10,'handlevisibility', 'off','MarkerEdgeColor','k','MarkerFaceColor',hmap(NoA+ns,:))];
end
plot(nan,nan,'p','MarkerSize',10,'handlevisibility', 'on','MarkerEdgeColor','k','MarkerFaceColor','none');
plot(ts.x,ts.y,'ko','MarkerSize',6,'MarkerFaceColor','k');
plot(Domain,'FaceColor','none','EdgeColor', 'b','HandleVisibility','off');
axis([xmin xmax ymin ymax]);
xlabel('$x$ (m)'); ylabel('$y$ (m)');
axis square
hold off;
legendlabelsS{end+1} = 'Est.';
legendlabelsS{end+1} = 'Source';
legend(legendlabelsS,'Orientation','horizontal','FontSize',12,'Location','bestoutside','Position',[0.440506315221754,0.614356939356939,0.425755380794701,0.028311965811966]);


timetxt = ['Time = ',num2str(tnow),'s and No. of Sources: ',num2str(length(EstHistory(end).estconfidence))];
title(at1,timetxt,'FontSize',12,'Interpreter','latex');
drawnow
end