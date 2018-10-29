function y = SoftThresh(x,t)
% y = SoftThresh(x,t)
% Soft-thresholding function. Works for both real or complex-valued valued 
% inputs 'x'
s = abs(x) - t;
s = (s + abs(s))/2;
y = sign(x).*s; %
    