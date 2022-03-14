function [FNP] = IQR_FNP(sig)
%UNTITLED7 Summary of this function goes here
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
sig_xx = sig(:,1);
sig_yy = sig(:,2);
tau_xy = sig(:,3);

%--------------------------------------------------------------------------
%                  Berechnung Nichtproportionalitätskennzahl
%--------------------------------------------------------------------------

% Hauptspannungsberechnung 

s2 = 0.5 * (sig_xx + sig_yy) - sqrt( 0.25*(sig_xx - sig_yy).^2 + tau_xy.^2 ); 

% Berechnung der Hauptspannungsrichtung

PSI = atan(tau_xy./(sig_xx - s2));
PSI(isnan(PSI)) = 0;

% Berechnung der Quantile

PSI_SORTED = sort(PSI);

n = size(PSI_SORTED,2);
p25 = 0.25;
p75 = 0.75;
if ( (n*p25) == floor(n*p25) )
    PSI25 = PSI_SORTED(n*p25);
    
else
    PSI25 = PSI_SORTED(floor(n*p25));
    
end

if ( (n*p75) == floor(n*p75) )
    PSI75 = PSI_SORTED(n*p75);
    
else
    PSI75 = PSI_SORTED(floor(n*p75));
    
end

IQR = PSI75 - PSI25;

FNP = IQR/(pi/2);
end

