function [FNP] = MeggiolaroFNP(sig)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
% Aus Masterarbeit Andreas Kern

%--------------------------------------------------------------------------
%                       Eingabedaten überprüfen
%--------------------------------------------------------------------------

if size(sig,2) ~= 3
    msg = 'kein ebener Spannungszustand';
    error(msg)
end

%--------------------------------------------------------------------------
%                       Spannungen identifizieren
%--------------------------------------------------------------------------
sig_xx = sig(:,1)';
sig_yy = sig(:,2)';
tau_xy = sig(:,3)';

%--------------------------------------------------------------------------
%                  Berechnung Nichtproportionalitätskennzahl
%--------------------------------------------------------------------------
X = sig_xx - sig_yy;
Y = sqrt(3)*tau_xy;

Xtk = X(1:end-1)';
Xtk1 = X(2:end)';
Xci = Xtk + (Xtk1 - Xtk)/2;

Ytk = Y(1:end-1)';
Ytk1 = Y(2:end)';
Yci = Ytk + (Ytk1 - Ytk)/2;

dpi = sqrt( (Xtk1 - Xtk).^2 + (Ytk1 - Ytk).^2 );
dpi(dpi == 0) = 1e-30;
psi = asin( (Ytk1 - Ytk)./ dpi );

Ixx = 1/sum(dpi) * sum( (dpi.^2/12 .* sin(psi).^2 + Yci.^2 ) .* dpi );
Iyy = 1/sum(dpi) * sum( (dpi.^2/12 .* cos(psi).^2 + Xci.^2 ) .* dpi );
Ixy = -1/sum(dpi) * sum( (dpi.^2/12 .* sin(psi).*cos(psi)+Yci.*Xci ) .* dpi );

l1 = (Ixx + Iyy)/2 + sqrt( (Ixx - Iyy)^2/4 + Ixy^2 );
l2 = (Ixx + Iyy)/2 - sqrt( (Ixx - Iyy)^2/4 + Ixy^2 );

FNP = sqrt( l2/l1 );
end

