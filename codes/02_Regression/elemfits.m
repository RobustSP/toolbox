function [beta, w] = elemfits(y,x) 
% [beta, w] = elemfits(x,y) 
% elemfits compute the Nx(N-1)/2 elemental fits, i.e., intercepts b_{0,ij}
% and slopes b_{1,ij}, that define a line y = b_0+b_1 x that passes through 
% the data points (x_i,y_i) and (x_j,y_j), i<j, where i, j in {1, ..., N}. 
% and the respective weights | x_i - x_j | 
% INPUT: 
%    y : (numeric) N x 1 vector of real-valued outputs (response vector)
%    x : (numeric) N x 1 vector of inputs (feature vector) 
% OUTPUT:
%  beta: (numeric) N*(N-1)/2 matrix of elemental fits 
%     w: (numeric) N*(N-1)/2 matrix of weights
%
% version: Aug 31, 2018 (Dependencies: uses ladlasso function
% authors: Esa Ollila 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N=length(x); 
B=repmat((1:N)',1,N);
A=B';
a=A(A<B); 
b=B(A<B);
w = (x(a) - x(b));
beta = zeros(2,length(a));
beta(2,:) = (y(a) - y(b))./w;  
beta(1,:) = (x(a).*y(b) - x(b).*y(a))./w;
w = abs(w); 

