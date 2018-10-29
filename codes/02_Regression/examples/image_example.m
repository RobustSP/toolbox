clear; clc;
addpath('toolbox');
addpath('toolbox/data/');

%- Read image of sqares (20 x 20 pixels) 
load images.mat % contains the vectors y20 (clean data) and y20n (noisy data)
n = numel(y20)
scaledata = @(x) 3*(x - min(x))./(max(x)-min(x));

%-- Plot the the image 
figure(1);clf; 
subplot(2,2,1)
imagesc(reshape(y20,sqrt(n),sqrt(n))); 
colormap gray;
axis square off; 
title('original image');

%-- Plot the signal 
figure(2); clf;
subplot(2,2,1);
plot(1:n,y20,'k*','MarkerSize',14); 
axis tight;
xlabel('$i$','Interpreter','Latex')
ylabel('$\beta_i$','FontSize',20,'Interpreter','Latex')
set(gca,'YTick',[0 1 2 3],'XTick',[1 100 200 300 400]);
set(gca','LineWidth',3,'FontSize',20,'FontName','TimesNewRoman');
title('original signal')

%-- Plot the image + noise 
figure(1)
subplot(2,2,2)
imagesc(reshape(y20n,sqrt(n),sqrt(n))); 
axis square off;
colormap gray;
title('image + noise');

%-- Plot the noisy (measured) signal 
figure(2); 
subplot(2,2,2);
plot(1:n,scaledata(y20n),'k*','MarkerSize',14); 
xlabel('$i$','FontSize',16,'Interpreter','Latex')
ylabel('$y_{i}$','FontSize',20,'Interpreter','Latex')
set(gca,'YTick',[0 1 2 3],'XTick',[1 100 200 300 400]);
set(gca,'LineWidth',3,'FontSize',20,'FontName','TimesNewRoman');
title('measured signal')


%-- Compute the Lasso solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 20; % grid size
[Blas20n,stats] = enetpath(y20n,eye(n),1,L,10^-3,false);
Blas20n = Blas20n(:,2:end); % get rid of all zeros (first column)
% Choose the best Lasso solution
ero = scaledata(Blas20n) - repmat(y20,1,size(Blas20n,2));
[MSElasso, indx] = min(sum(ero.^2)); MSElasso
Blas = Blas20n(:,indx); % the best Lasso solution 
lam_las = stats.Lambda(1+indx); % the best lambda value

figure(1); 
subplot(2,2,3);
imagesc(reshape(Blas,20,20)); 
axis square off; 
colormap gray; 
title(['Lasso: lambda = ', num2str(lam_las)]);

figure(2);
subplot(2,2,3);
plot(1:n,scaledata(Blas(:)),'k*','MarkerSize',14); 
axis tight;
xlabel('$i$','FontSize',16,'Interpreter','Latex')
ylabel('$\hat \beta_{i}$','FontSize',20,'Interpreter','Latex')
set(gca,'YTick',[0 1 2 3],'XTick',[1 100 200 300 400]);
set(gca,'LineWidth',3,'FontSize',20,'FontName','TimesNewRoman');
title('Lasso');

%-- Compute the Rank-FLasso solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start with some initial values of lambda1 and lambda2
lambda2 = 340;
lambda1 = 124; 
B1 = rankflasso(y20n,eye(n),lambda1,lambda2,Blas,1);
MSE_rank1 = sum((scaledata(B1)- y20).^2)

%-- adjust the parameters 
lambda2 = 420;
lambda1 = 35; 
B2= rankflasso(double(y20n),eye(20*20),lambda1,lambda2,B1,1);
MSE_rank2 = sum((scaledata(B2)- y20).^2)

figure(1); 
subplot(2,2,4)
imagesc(reshape(B2,20,20)); 
axis off square; 
title(['Rank-FLasso: lambda1= ', num2str(lambda1), ' lambda2 = ', num2str(lambda2)]);
colormap gray; 

figure(2);
subplot(2,2,4);
plot(1:n,scaledata(B2),'k*','MarkerSize',14); 
axis tight;
xlabel('$i$','FontSize',16,'Interpreter','Latex')
ylabel('$\hat \beta_{i}$','FontSize',20,'Interpreter','Latex')
set(gca,'YTick',[0 1 2 3],'XTick',[1 100 200 300 400]);
set(gca,'LineWidth',3,'FontSize',20,'FontName','TimesNewRoman');
title('Rank-FLasso');