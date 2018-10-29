function y=eta(x,c)
 if nargin<2
    c = 1; 
 end
x = x/c;
 
y(abs(x)<=2)=x((abs(x)<=2));
y(and(abs(x)>2,abs(x)<=3))=0.016*x(and(abs(x)>2,abs(x)<=3)).^7-0.312*x(and(abs(x)>2,abs(x)<=3)).^5+1.728*x(and(abs(x)>2,abs(x)<=3)).^3-1.944*x(and(abs(x)>2,abs(x)<=3));
y(abs(x)>3)=0;
y=c*y';
end
