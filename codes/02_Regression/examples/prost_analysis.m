clear; clc;
addpath('/examples');

load data/prostate
[n, p] = size(X);
Xone = [ones(n,1) X];
LSE = Xone \ y % Least squares estimate
GRlen = 120;

%-- LASSO 
%%%%%%%%%%%
[B,stats] = enetpath(y,X,1);
[~,k] = min(stats.BIC); 
blas = B(:,k); % LASSO-BIC solution
bmaxlas = sum(abs(B(2:end,end))); % largest value of || \beta ||_1
%- PLOT 
h=figure(1); clf;
subplot(2,2,1);
locs = B(2:end,end);
locs(3) = locs(3)-0.025; % 'age' is too close, so put it down
locs(7) = locs(7)+0.01;  % 'gleason' is too close, so put it up
loc_x = sum(abs(blas(2:end)))/bmaxlas;
prostate_plot_setup(sum(abs(B(2:end,:)))/bmaxlas,B(2:end,:),locs,loc_x,names)

%-- LAD-LASSO 
%%%%%%%%%%%%%%
reltol = 1.0e-7;
tic;[Blad, statslad] = ladlassopath(y, X,[],[],[],reltol); 
[~,ladind] = min(statslad.gBIC);
blad = Blad(:,ladind); % LAD-Lasso BIC solution
bmaxlad = max(sum(abs(Blad(2:end,:)))); % largest solution || \beta ||_1

%- PLOT
figure(1); subplot(2,2,3); 
locs = Blad(2:end,end);
locs(2) = locs(2)+0.02; % lweight up
locs(7) = locs(7)+0.02; % gleason up
locs(3) = locs(3)-0.02; % age down
loc_x = sum(abs(blad(2:end)))/bmaxlad;
prostate_plot_setup(sum(abs(Blad(2:end,:)))/bmaxlad,Blad(2:end,:),locs,loc_x,names)

%-- Rank-LASSO
%%%%%%%%%%%%%%%
[Brlad,~,statsrlad] = ranklassopath(y, X); 
[~,rladind] = min(statsrlad.gBIC);
brlad = Brlad(:,rladind); 
bmaxrlad = max(sum(abs(Brlad)));

%- PLOT
figure(1); 
subplot(2,2,4)
locs = Brlad(:,end);
locs(7) = locs(7)+0.01; % gleason up
locs(3) = locs(3)-0.03; % age down
loc_x = sum(abs(brlad))/bmaxrlad;
prostate_plot_setup(sum(abs(Brlad))/bmaxrlad,Brlad,locs,loc_x,names)

%-- M-LASSO
%%%%%%%%%%%%
[Bhub,~,statshub] = hublassopath(y,X);
[~,hubind] = min(statshub.gBIC);
bhub = Bhub(:,hubind); 
bmaxhub = max(sum(abs(Bhub)));

%- PLOT
figure(1); subplot(2,2,2); 
locs = Bhub(:,end);
loc_x = sum(abs(bhub))/bmaxhub;
locs(3) = locs(3)-0.03; % age down
locs(7) = locs(7)+0.02; % gleason up
prostate_plot_setup(sum(abs(Bhub))/bmaxhub,Bhub,locs,loc_x,names)

%-- OUTLIER --       
%%%%%%%%%%%%%%

yout = y;
yout(1) = 55;

%-- LASSO 
%%%%%%%%%%%
[Bout,stats2] = enetpath(yout,X,1);
[~,k] = min(stats2.BIC); 
blas_out = Bout(:,k); % LASSO-BIC solution
bmaxlas_out = sum(abs(Bout(2:end,end))); % largest value of || \beta ||_1

%- PLOT 
figure(2); clf;
subplot(2,2,1);
locs = Bout(2:end,end);
locs(1) = locs(1)-0.04; % lcacvol down
locs(3) = locs(3)-0.08; % age down
locs(4) = locs(4)+0.06; % lbph up
locs(6) = locs(6)+0.02; % lcp up
locs(8) = locs(8)+0.07; % pgg45 up
loc_x = sum(abs(blas_out(2:end)))/bmaxlas_out;
prostate_plot_setup(sum(abs(Bout(2:end,:)))/bmaxlas_out,Bout(2:end,:),locs,loc_x,names,[])
    
%-- LAD-LASSO 
%%%%%%%%%%%%%%
reltol = 1.0e-7;
tic;[Blad2, statslad2] = ladlassopath(yout,X,[],[],[],reltol); 
[~,ladind2] = min(statslad2.gBIC);
blad2 = Blad2(:,ladind2); % LAD-Lasso BIC solution
bmaxlad2 = max(sum(abs(Blad2(2:end,:)))); % largest solution || \beta ||_1

%- PLOT
figure(2); subplot(2,2,3); 
locs = Blad2(2:end,end);
locs(6) = locs(6)-0.04; % lcp 
locs(8) = locs(8)+0.02; % pgg45 up
locs(7) = locs(7)+0.02; % gleason up
locs(3) = locs(3)-0.02; % age down
loc_x = sum(abs(blad2(2:end)))/bmaxlad2;
prostate_plot_setup(sum(abs(Blad2(2:end,:)))/bmaxlad2,Blad2(2:end,:),locs,loc_x,names)

%-- Rank-LASSO
%%%%%%%%%%%%%%%
[Brlad2,~,statsrlad2] = ranklassopath(yout,X); 
[~,rladind2] = min(statsrlad2.gBIC);
brlad2 = Brlad2(:,rladind2); 
bmaxrlad2 = max(sum(abs(Brlad2)));

%- PLOT
figure(2); 
subplot(2,2,4)
locs = Brlad2(:,end);
locs(7) = locs(7)+0.04; % gleason up
locs(3) = locs(3)-0.025; % age down
locs(6) = locs(6)-0.015; % lcp 
loc_x = sum(abs(brlad2))/bmaxrlad2;
prostate_plot_setup(sum(abs(Brlad2))/bmaxrlad2,Brlad2,locs,loc_x,names)

%-- M-LASSO
%%%%%%%%%%%%
[Bhub2,~,statshub2] = hublassopath(yout,X);
[~,hubind2] = min(statshub2.gBIC);
bhub2 = Bhub2(:,hubind2); 
bmaxhub2 = max(sum(abs(Bhub2)));

%- PLOT
figure(2); subplot(2,2,2); 
loc_x = sum(abs(bhub2))/bmaxhub2;
locs = Bhub2(:,end);
locs(3) = locs(3)-0.03; % age down
locs(7) = locs(7)+0.02; % gleason up
prostate_plot_setup(sum(abs(Bhub2))/bmaxhub2,Bhub2,locs,loc_x,names)