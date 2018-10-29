function psix = psihub(x,c)
% Huber's score function psi: input is N x 1 data vector x which can be complex
% or real and threshold contant c
psix = ( x.*(abs(x)<=c) + c*sign(x).*(abs(x)>c));
