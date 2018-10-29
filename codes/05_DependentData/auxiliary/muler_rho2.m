function rho=muler_rho2(x)    

rho(abs(x)<=2)=0.5*x((abs(x)<=2)).^2;
rho(and(abs(x)>2,abs(x)<=3))=0.002*x(and(abs(x)>2,abs(x)<=3)).^8-0.052*x(and(abs(x)>2,abs(x)<=3)).^6+0.432*x(and(abs(x)>2,abs(x)<=3)).^4-0.972*x(and(abs(x)>2,abs(x)<=3)).^2+1.792;
rho(abs(x)>3)=3.25;
rho=rho';
end
