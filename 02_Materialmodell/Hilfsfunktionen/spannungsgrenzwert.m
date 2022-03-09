function siggrenz = spannungsgrenzwert(material,para,M)
% Funktion berechnet Spannungsgrenzwert der einzelnen Modelle
% INPUT:
%  material - str mit modellname
%  M        - ANzahl Backstresstensoren
%  para     - Materialparameter
% OUTPUT:
% siggrenz  - (einachsiger) Spannungsgrenzwert
% -------------------------------------------------------------------------

% Radius FF
r0 = para(end);

% unterscheide Materialmodelle
switch material
    case 'Chaboche'
        ri = para(5+M:end-1);
        siggrenz = r0 +  sum(ri);
    case 'KarimOhno'
        ri = para(3+M:2+2*M);
        siggrenz = r0 +  sum(ri);
    case 'OhnoWang'
        ri = para(3+M:2+2*M);
        siggrenz = r0 + sqrt(3/2)*sum(ri);
    case 'Jiang'
        ri = para(9+5*M:9+6*M-1);
        siggrenz = r0 + sqrt(3/2)*sum(ri);
    case 'OWT'
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!
        % Nichtproportionale Verfestigung nicht berücksichtigt
        ri = para(3+M:2+2*M);
        siggrenz = r0 + sqrt(3/2)*sum(ri);
    otherwise
        msg = 'Modell nicht implementiert';
        error(msg);
end