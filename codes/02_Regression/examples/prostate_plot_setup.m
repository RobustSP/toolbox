function prostate_plot_setup(x,Y,locs,loc_x,names,axval)

hold on;
tmp = 0:0.1:0.7;
for ii =1:8
    plot(x,Y(ii,:),'Color',repmat(tmp(ii),1,3),'LineWidth',2.0)
end
axgca = gca;
axgca.GridLineStyle = '--';
if nargin < 6 
    axis([0 1.2 -0.25 0.9])
    set(gca,'YTick',[-0.2 0 0.2 0.4 0.6 0.8 1.0],'FontSize',12)
end
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1.0],'FontSize',12)

text(repmat(1.02,8,1),locs,names(1:8),'FontSize',14);
plot([loc_x loc_x],[-0.25 1.05],'k--','LineWidth',2.0)

ylabel('Coefficients');
xlabel('normalized $\|\hat \beta(\lambda) \|_1$','Interpreter','Latex')
grid on
set(gca,'LineWidth',3,'FontName','Helvetica','FontSize',18);