function rhox = rhotuk(x,c)
rhox = (c^2/3)*((1-(1-(abs(x)/c).^2).^3).*(abs(x)<=c) + (abs(x)>c) ); % Tukey 
end


