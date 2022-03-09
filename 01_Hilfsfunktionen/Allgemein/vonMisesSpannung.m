function SIGV = vonMisesSpannung(SIG)
% Funktion berechnet von Mises Vergleichsspannung 
%
% INPUT:
% SIG      - Spannungstensor (ntens x ndata)
%
% OUTPUT:
% SIGV     - Vergleichsspannung (1 x ndata
% -------------------------------------------------------------------------

% Abbildungen
ntens = size(SIG,1);
if ntens == 6
    ndi = 3;
elseif ntens == 3
    ndi = 2;
elseif ntens == 2
    ndi = 1;
else
    msg = 'Falscher Spannungszustand';
    error(msg)
end
[P,PLINE] = set_maps(ntens,ndi);

% Spannungsdeviator
if ntens == 2
    S = P .* SIG;
else
    S = P * SIG;
end

% Vergleichsspannung
if ntens == 2
    SIGV = sqrt(3/2) * sqrt(sum( (PLINE.*S) .* S));
else
    SIGV = sqrt(3/2) * sqrt(sum( (PLINE*S) .* S));
end


end