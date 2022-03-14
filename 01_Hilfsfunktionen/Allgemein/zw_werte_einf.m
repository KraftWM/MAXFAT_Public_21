function [ukp_mod] = zw_werte_einf(nzw, ukp)
% Funktion fuegt linear Zwischenpunkte zu gegebener Umkehrpunktfolge hinzu
% INPUT: 
% nzw - Anzahl Zwischenwerte
% ukp - Umkehrlastpunktfolge
% OUTPUT:
% ukp_mod - modifizierte Umkehrlastpunktfolge
%__________________________________________________________________________

% Anzahl Datenpunkte
ndp = size(ukp,2);

% Abzahl Lastkanäle
nk = size(ukp,1);

ukp_mod = zeros(nk,ndp*(nzw+1)-nzw);
zwp = 0 : 1/(nzw+1) : 1-1/(nzw+1);
eins = ones(1,nzw+1);
zeiger = 1;
for i = 1:ndp-1
    d = ukp(:,i+1)-ukp(:,i);
    zw = ukp(:,i).*eins + d .* zwp;
    ukp_mod(:,zeiger:zeiger+nzw) = zw;
    zeiger = zeiger + nzw + 1;
end
ukp_mod(:,end) = ukp(:,end);