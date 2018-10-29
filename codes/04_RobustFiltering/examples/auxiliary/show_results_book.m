
fig1 = figure(1);
fig2 = figure(2);
fig3 = figure(3);      


% plotting trajectories
% Least-squares estimation
for ii = 1:parameter.N
ekf_th_x(ii) = ekf_th(1,1,ii);
ekf_th_y(ii) = ekf_th(1,2,ii);
end



% Robust M-estimation
for ii = 1:parameter.N
ekf_Hc_x(ii) = ekf_Hc(1,1,ii);
ekf_Hc_y(ii) = ekf_Hc(1,2,ii);
end

parameter


figure(1)
plot(parameter.BS(:,1),parameter.BS(:,2),'ok','Linewidth',2)
hold on
plot(parameter.thx(end,:),parameter.thy(end,:),'Linewidth',2)   
plot(ekf_th_x,ekf_th_y,'Linewidth',2)
plot(ekf_Hc_x,ekf_Hc_y,'Linewidth',2)

% evaluation of MSE

figure(2)
indices =sum(parameter.noiseindices)';
stem(find(indices>0),800*indices(indices>0),'y', 'Color',[ 0.9297    0.8016     0.5],'Linewidth',2,'Markersize',1) 
hold on
[eval_ekf ] = eval_track(ekf_th,parameter,'r',fig1,fig2,fig3);
[eval_ekf_Hc ] = eval_track(ekf_Hc,parameter,[.2 .1 .5],fig1,fig2,fig3);
        
        
        figure(3)
        legend('EKF','Robust EKF')
        
        figure(2)
        legend('NLOS','EKF','Robust EKF')    
        axis([0 parameter.N 0 800])
        
        figure(1)
        legend('BS','true','EKF','Robust EKF')  


