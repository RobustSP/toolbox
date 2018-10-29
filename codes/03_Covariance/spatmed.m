function  smed = spatmed(X,printitn)
%  Computes the spatial median based on (real or complex) data matrix X.
%  INPUT:
%         X: Numeric data matrix of size N x p. Each row represents one 
%           observation, and each column represents one variable 
% printitn : print iteration number (default = 0, no printing)
%
% OUTPUT
%      smed: Spatial median estimate 
%
% version: Sep 2, 2018 
% authors: Esa Ollila 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

len = sum(X.*conj(X),2); 
X = X(len~=0,:);
n = size(X,1);

if nargin <  2
    printitn = 0;
end

if isreal(X)
    smed0 = median(X);
else
    smed0 = mean(X);
end
norm0 = norm(smed0);

iterMAX = 500;
EPS = 1.0e-6;
TOL = 1.0e-5;   

for iter = 1:iterMAX 

   Xc = bsxfun(@minus,X,smed0);
   len = sqrt(sum(Xc.*conj(Xc),2)); 
   len(len<EPS)= EPS;
   Xpsi = bsxfun(@rdivide, Xc, len);
   update = sum(Xpsi)/sum(1./len);
   smed = smed0 + update; 
   
   dis = norm(update)/norm0;  
  
   if mod(iter,printitn)==0
        fprintf('At iter = %3d, dis=%.7f\n',iter,dis);
   end
   
   if (dis<=TOL) 
       break;             
   end
   smed0 = smed;
   norm0 = norm(smed);
   
end
% DONE