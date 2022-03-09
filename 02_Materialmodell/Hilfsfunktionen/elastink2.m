function xel = elastink2(s,a,P2,r,ds,Ftol,varargin)
% Funktion zum ermitteln des elastischen Anteils einer Belastung für eine
% Fließfunktion nach von Mises F = (s-a):(s-a)-2/3r^2 = 0
% INPUT:
%       s    -> Spannungsdeviator
%       a    -> Backstresstensor
%       P2   -> Skalarprodukt Spannungsdeviatoren
%       r    -> Radius FF
%       ds   -> Deviatorisches Spannungsinkrement
%   n_normal -> Anzahl normalenkomponenten
%   Ftol     -> Toleranz für Fließfläche
%   varargin -> Variabler Input (Bsp. Jiang Modell hat andere FF)
%              'jiang'  -> jiang modell
%              'doring' -> döring modell
%
%
% Tensoren werden als Vektoren der Dimension R^(ntens x 1) eingegeben
% Normalenkomponenten des Tensors sind an den ersten n_normal stellen des
% vektors
% OUTPUT:
%       xel -> elastischer Anteil des Inkrements
%              in [0,1]
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
% |  Stand: Juni 2019                                                    |
%  ----------------------------------------------------------------------


% Überprüfe Input
ntens = size(s,1);

% Abfangen von Nullen
sums = sum(abs(ds) <= 1e-10);
if sums == ntens
    xel = 1;
    return
end

% Sonderfall Jiang/Döring Modell
% F = (s-a):(s-a)-2/3r^2 = 0       -> r = einachs. FließNORMALspann.
% F = (s-a):(s-a)-2k^2 = 0 (jiang) -> k = einachs. FließSCHUBspann.
% 2 k^2 = 2/3 r^2 -> r = sqrt(3) k
% F = (s-a):(s-a)- k^2 = 0 (döring) -> k = deviatorische einachs. Fließspan.
% r = sqrt(3/2) k
if nargin > 6
    if strcmp(varargin{1,1},'jiang')
        r = sqrt(3) * r;
    elseif strcmp(varargin{1,1},'doring')
        r = sqrt(3/2) * r;
    end
end

% Relativspannung
b = s - a; 

% Unterscheide Spannungsszutände
if ntens == 1        % 1D
    nds2 = ds^2;                         % Normquadrat Dev.Spannungsinkrement
    q1 = b * ds;                         % quasi der Belastungsfaktor
    q2 = b^2-r^2;                        % Aktuelle Fließfläche
elseif ntens == 2    % Sigma - Tau
    nds2 = ds' * (P2 .* ds);
    q1 = b' * (P2 .* ds);
    q2 = b' * (P2 .* b) - 2/3 * r^2;
else                 % 3D & ESZ
    nds2 = ds' * P2 * ds;
    q1 = b' * P2 * ds;
    q2 = b' * P2 * b - 2/3 * r^2;
end

% Begrenze Fließspannung
if (q2>0 && q2<Ftol) || (q2<0 && q2>-Ftol)
        q2 = 0;
end


% Berechne elastioschen Anteil
x1 = (-q1 + sqrt(q1^2 - q2*nds2))/nds2;
x2 = (-q1 - sqrt(q1^2 - q2*nds2))/nds2;
xel = max(x1,x2);

% Abbruch bei Fehler
xtol = 1e-4; % setze Toleranz
if xel < 0
    % kleiner fehler wird korrigiert
    if xel > -xtol 
        xel = 0;
    % Bei großem Fehler Rechnung abbrechen    
    else
        % Mögliche Ursache ist das Konsistenzbedingung nicht eingehalten wurde
        msg='elastischer Anteil kleiner Null berechnet';
        error(msg)
    end
elseif abs(imag(xel)) > 0
    msg=['imaginärer elastischer Teilschritt, ',...
         'Ursache ist zu große Überspannung ', 'q1 = ', num2str(q1), ...
         ';  q2 = ' , num2str(q2), '. q2 ist immer kleiner gleich Null',...
         ' solange keine Überspannung auftritt'];
    error(msg)
elseif xel > 1 
    % kleine Fehler Korrigieren
    if xel < xel + xtol
        xel = 1;
    % Bei großem Fehler Rechnung abbrechen    
    else
        msg='elastischer Anteil größer 1 berechnet';
        error(msg)
    end
end
 