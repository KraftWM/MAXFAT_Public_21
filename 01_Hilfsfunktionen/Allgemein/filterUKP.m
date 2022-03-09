% -------------------------------------------------------------------------
% Hilfsfunktion filtern Umkehrpunkte
function [UKP,IDX] = filterUKP(L)
% Bilde 1. Ableitung
dL = diff(L,[],2);
% Vorzeichenwechsel
UKP = dL(:,1:end-1).*dL(:,2:end) <= 0 & L(:,3:end) ~= L(:,2:end-1);
UKP = [true,UKP,true]; % Ersten und letzten Punkt erstmal behalten
% Prüfe 1. Punkt !! Annahme Lastfolge startet in 0
if sign(L(2)) == sign(L(1))         % Nur Filtern wenn Gleiches VZ
    if abs(L(2)) > abs(L(1))        % Wenn 1. Punk näher am Ursprung liegt
        UKP(1) = false;
    end
end
% Indices der UKP
IDX = 1:size(L,2);
IDX = IDX(UKP)';
end % Ende UKP Filter