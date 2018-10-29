function rhox = rhohub(x,c)
% Huber's loss function rho: input is N x 1 data vector x which can be complex
% or real and threshold contant c
rhox =  (abs(x).^2).*(abs(x)<=c) + (2*c*abs(x)-c^2).*(abs(x)>c);
end

