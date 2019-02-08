function wx = whub(absx,c)
% Huber's weight function w: input is N x 1 data vector x which can be complex
% or real and threshold contant c
wx = 1.*(absx<=c) + (c*(1./absx)).*(absx>c);
wx(absx==0) = 1;
end
