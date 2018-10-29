clear;
NR_ITER = 10000; 
p = 8;  % dimension 
s = exp(1i*(1:p)*pi).'; %  pulse of || p ||^2 = m
alpha = [0.001:0.002:0.04 0.04]; % set of PFA's 

lam =  1- alpha.^(1/(p-1));  % = betainv(1-alpha,1,m-1); 
nu = 0.5;
nlist = [4*p 8*p 12*p];  % number of secondary dta
ntest = 10; % number of primary data 

Pfa1 = zeros(numel(nlist),numel(alpha));
Pfa2 = zeros(numel(nlist),numel(alpha));

Lambda1 = zeros(NR_ITER,ntest);
Lambda2 = zeros(NR_ITER,ntest);

rng default % For reproducibility

for kk = 1:numel(nlist)
    
    n = nlist(kk);

    for it = 1:NR_ITER 

        % Generate covariance matrix
        %l = rand(1,p);
        l = unifrnd(0.1,1,[1 p]);
        P = orth(rand(p,p) + 1i*rand(p,p));
        sig  = p*P*diag(l./sum(l))*P';
        sig(1:(p+1):p^2)=real(sig(1:(p+1):p^2)); 
        sqrsig = sqrtm(sig);
    
        % Generate the secondary data and compute the covariance 
        x0 = sqrsig*sqrt(1/2)*(randn(p,n) + 1i*randn(p,n)); % ~ N_p(0,I)
        x = repmat(sqrt(gamrnd(nu,1/nu,1,n)),p,1).*x0 ; % ~ K_p,v(0,I)
        hsig1 = Mscat(x','t-loss',0);   % compute Tyler's M-estimator
        hsig2 = Mscat(x','Huber',0.8);  % Huber's M-estimator
        B1 = sqrtm(inv(hsig1)); % 
        B2 = sqrtm(inv(hsig2)); % 

        % Generate primary data from C K_v(0,sig);
        z0 = sqrsig*sqrt(1/2)*(randn(p,ntest)+ 1i*randn(p,ntest)); 
        z = repmat(sqrt(gamrnd(nu,1/nu,1,ntest)),p,1).*z0 ;% primary data   
        
        % Compute the ADAPTIVE DETECTOR 
        v1 = B1*z; q1 = B1*s; q1 = q1/norm(q1); 
        v2 = B2*z; q2 = B2*s; q2 = q2/norm(q2); 
        v1 = bsxfun(@rdivide, v1, sqrt(sum(v1.*conj(v1))));
        v2 = bsxfun(@rdivide, v2, sqrt(sum(v2.*conj(v2))));       
        Lambda1(it,:) = abs(v1'*q1).^2; Lambda2(it,:) = abs(v2'*q2).^2;

        if mod(it,250) ==0, fprintf('. '); end
    end
    fprintf('\n');
    for al = 1:numel(alpha) 
       Pfa1(kk,al) = sum(Lambda1(:) > lam(al))/(NR_ITER*ntest);
       Pfa2(kk,al) = sum(Lambda2(:) > lam(al))/(NR_ITER*ntest);
    end
end

beep;

%-- Plot  Threshold  vs PFA
figure(2); clf
subplot(2,2,1);
plot(lam,alpha,'LineStyle','-','Color','k','LineWidth',2)
axis([min(lam) max(lam) 0 0.08]);  hold on; 
plot(lam,Pfa1(3,:),'--',lam,Pfa1(2,:),':',lam,Pfa1(1,:),'-.','LineWidth',2)
set(gca,'YTick',[0 0.01 0.03 0.05 0.07])
ylabel('$P_{\mathrm{FA}}$','Interpreter','Latex')
xlabel('Detection threshold $(\lambda)$','Interpreter','Latex')
hleg=legend('$n=\infty$','$n=96$','$n=64$','$n=32$');
set(hleg,'Interpreter','Latex','FontSize',20)
grid on
set(gca,'LineWidth',3,'FontSize',18);

subplot(2,2,2);
plot(lam,alpha,'LineStyle','-','Color','k','LineWidth',2)
axis([min(lam) max(lam) 0 0.08]);  hold on; 
plot(lam,Pfa2(3,:),'--',lam,Pfa2(2,:),':',lam,Pfa2(1,:),'-.','LineWidth',2)
set(gca,'YTick',[0 0.01 0.03 0.05 0.07])
ylabel('$P_{\mathrm{FA}}$','Interpreter','Latex')
xlabel('Detection treshold $(\lambda)$','Interpreter','Latex')
hleg=legend('$n=\infty$','$n=96$','$n=64$','$n=32$');
set(hleg,'Interpreter','Latex','FontSize',20)
grid on
set(gca,'LineWidth',3,'FontSize',18);

%-- Plot the histogram 
subplot(2,2,3); 
histogram(Lambda1,30,'normalization','pdf')
ylabel('Frequency','FontSize',20,'FontName','Helvetica');
xlabel('$\hat \Lambda$','FontSize',20,'FontName','Helvetica','Interpreter','Latex');
hold on 
xval = 0:0.01:1;
plot(xval,betapdf(xval,1,p-1),'k-','LineWidth',2)
axis([0 1 0 max(betapdf(xval,1,p-1))])
set(gca,'LineWidth',3,'FontSize',18);
grid on

%--  Q-Q plot for Tyler's M-estimator for n = 96
lams = sort(Lambda1(:));  % samples
len = numel(lams);
pvals = ((1:len) - (1/2))/len;
qvals = betainv(pvals,1, p-1); % quantiles of exp with mu=1
q99 = betainv(0.99,1,p-1);
q999 = betainv(0.999,1,p-1);
e99 = lams(len-1000); e999 = lams(len-100);

subplot(2,2,4);
plot(qvals,lams,'.','Color',repmat(0.2,1,3),'LineWidth',1.4); 
hold on;
plot([0 1],[0 1],'k--','LineWidth',1.4) % straight line
axis([0 0.85 0 0.85])
set(gca,'YTick',0:0.1:0.8,'Xtick',0:0.1:0.8)

ylabel('Empirical quantile','FontSize',20,'FontName','Helvetica');
xlabel('Theorerical quantile','FontSize',20,'FontName','Helvetica');
plot([q99 q99], [0 q99],'k-'); plot([0 e99], [e99 e99],'k-')
plot([q999 q999], [0 q999],'k-'); plot([0 e999], [e999 e999],'k-')
text(0.05,0.55,'1%','FontSize',16); text(0.05,0.7,'0.1%','FontSize',16);
text(0.4,0.05,'1%','FontSize',16); text(0.52,0.05,'0.1%','FontSize',16);
set(gca,'LineWidth',3,'FontSize',18);