function [fnpBOL] = BolchounFNP(sig)

sig_xx = sig(:,1);
sig_yy = sig(:,2);
tau_xy = sig(:,3);

nT = 150;
theta = linspace(0,pi,nT);
dT = pi/nT;

coeff = zeros(1,nT);
for i = 1 : nT
    phi = theta(i);
    sxx = 0.5*(sig_xx + sig_yy) + (sig_xx - sig_yy)*cos(2*phi) + tau_xy*sin(2*phi);
    txy = 0.5*(sig_yy - sig_xx)*sin(2*phi) + tau_xy*cos(2*phi); 
    
    n = size(sxx,1);
    XY = sum(sxx .* txy,1);
    XS = sum(sxx,1);
    YS = sum(txy,1);
    XS2 = sum(sxx.^2,1);
    YS2 = sum(txy.^2,1);
    cor = (n*XY - XS*YS)/( sqrt((n*XS2 - (XS)^2)*(n*YS2 - (YS)^2)) );
    coeff(:,i) = cor;
    
end

M = 1/pi * sum(coeff(2:end).^2*dT,2);

fnpBOL =  1- M;

end


