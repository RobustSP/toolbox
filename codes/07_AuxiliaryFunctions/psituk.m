function psix = psituk(x,c)
psix = x.*((1-(abs(x)/c).^2).^2).*(abs(x)<=c);
end 
