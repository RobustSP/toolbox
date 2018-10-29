function wx = wtuk(absx,c)
wx = ((1-(absx/c).^2).^2).*(absx<=c);
end 
