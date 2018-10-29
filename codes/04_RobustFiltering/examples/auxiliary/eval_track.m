function [eval_filter rmse ] = eval_track(xx,param,color,fig1,fig2,fig3)
% evaluates the rmse, med, etc. of several tracks


% input variables:

% xx            % estimated state vector (mc x N)

% param         structure that contains the different tracks and the true
                % trajectory
                
% color         % indicates the color of the track
  
% fig           % figure handle of the figure the curves are plotted



% output variables

% eval_filter      structure that contains the rmse, med and other metrices


% if (nargin==5)
%     fig3 = figure(3);
% end


    % intitialization
    eex(param.mc,param.N)  = 0;
    eey(param.mc,param.N)  = 0;
    med(param.mc,param.N)  = 0;
    eevx(param.mc,param.N) = 0;
    eevy(param.mc,param.N) = 0;
bx(param.mc,param.N) = 0;
by(param.mc,param.N) = 0;
bvx(param.mc,param.N) = 0;
bvy(param.mc,param.N) = 0;

    
    
 for ii = 1: param.mc
    for kk=1:param.N
        bx(ii,kk)   =  (xx(ii,1,kk)-param.thx(ii,kk));
        by(ii,kk)   =  (xx(ii,2,kk)-param.thy(ii,kk));
        bvx(ii,kk)   =  (xx(ii,3,kk)-param.thvx(ii,kk));
        bvy(ii,kk)   =  (xx(ii,4,kk)-param.thvy(ii,kk));
        eex(ii,kk)  = (xx(ii,1,kk)-param.thx(ii,kk)).^2;
        eey(ii,kk)  = (xx(ii,2,kk)-param.thy(ii,kk)).^2;
        med(ii,kk)  = sqrt((xx(ii,1,kk)-param.thx(ii,kk)).^2 +(xx(ii,2,kk)-param.thy(ii,kk)).^2);
        eevx(ii,kk)  = (xx(ii,3,kk)-param.thvx(ii,kk)).^2;
        eevy(ii,kk)  = (xx(ii,4,kk)-param.thvy(ii,kk)).^2;
    end
 end
    
 
eval_filter.rmsex = sqrt(mean(mean(eex(:,param.discardN:end))));
eval_filter.rmsey = sqrt(mean(mean(eey(:,param.discardN:end))));

eval_filter.rmsevx = sqrt(mean(mean(eevx(:,param.discardN:end))));
eval_filter.rmsevy = sqrt(mean(mean(eevy(:,param.discardN:end))));


eexy = eex+eey;

eval_filter.med  = mean(mean(med(:,param.discardN:end)));
eval_filter.rmse = sqrt(mean(mean(eexy(:,param.discardN:end))));
eval_filter.biasx_v = mean(bx(:,param.discardN:end));
eval_filter.biasy_v = mean(by(:,param.discardN:end));
eval_filter.biasvx_v = mean(bvx(:,param.discardN:end));
eval_filter.biasvy_v = mean(bvy(:,param.discardN:end));
eval_filter.biasx = mean(mean(bx(:,param.discardN:end)));
eval_filter.biasy = mean(mean(by(:,param.discardN:end)));
eval_filter.biasvx = mean(mean(bvx(:,param.discardN:end)));
eval_filter.biasvy = mean(mean(bvy(:,param.discardN:end)));

 

MED = med(:,param.discardN:end);
MED1=med;
    
    thx(param.N) = 0;    
    thy(param.N) = 0;
    
    for kk=1:param.N
        thx(kk) = xx(end,1,kk);
        thy(kk) = xx(end,2,kk);
    end
        
    
    linewidth = 2;
    markersize = 12;
    
    if (param.figure)
        figure(fig1)
        hold on
        grid on
        plot(thx,thy,'Color',color,'Linewidth',linewidth,'MarkerSize',markersize)
        xlabel('x-position in m')
        ylabel('y-position in m')
        set(get(gcf,'CurrentAxes'),'FontSize',16);
        figure(fig2)
        
        if (length(mean(eexy))>1)
             rmse = sqrt(mean(eexy));
             kgrid = 1:param.grid:length(rmse);
         if (strcmp(param.plot,'mse'))
            %rmse = sqrt(mean(eexy));
            rmse_grid = rmse(1:param.grid:end);
            plot(kgrid,rmse_grid,'Color',color,'Linewidth',linewidth,'MarkerSize',markersize),
         elseif (strcmp(param.plot,'MED') )
             
             
             
            plot(kgrid,mean(MED1),'Color',color,'Linewidth',linewidth,'MarkerSize',markersize),
               
         end
        end
        if (length(mean(eexy))==1)
            plot(sqrt(eexy),'Color',color,'Linewidth',linewidth,'MarkerSize',markersize)
       % keyboard   
        
        end
        xlabel('Time')
        ylabel('Root Mean Square Error (RMSE)')
        %rmse
        set(get(gcf,'CurrentAxes'),'FontSize',16);
        hold on
        grid on
        figure(fig3)
        hold on
        [cdfy, cdfx ] = cdfcalc(MED(:));
        cdfy(end) = [];
        cdfy_grid = zeros(1,length(cdfy));
         
        cdfy_grid = cdfy(1:param.grid*param.mc:end);
        cdfx_grid = cdfx(1:param.grid*param.mc:end);
        
        set(get(gcf,'CurrentAxes'),'FontSize',16);
        ylabel('Empirical CDF','FontSize',16);
        xlabel('Error distance in m','FontSize',16)
       % keyboard
        h=plot(cdfx_grid,cdfy_grid,'Color',color,'Linewidth',linewidth,'MarkerSize',markersize);
     
        %   set(h,'fontsize',14)

        %cdfplot(MED(:))
        grid on        
    end




return



