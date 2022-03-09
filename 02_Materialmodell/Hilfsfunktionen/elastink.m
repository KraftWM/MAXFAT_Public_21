function xel = elastink(s,a,r,ds,n_normal,varargin)
% Funktion zum ermitteln des elastischen Anteils einer Belastung für eine
% Fließfunktion nach von Mises F = (s-a):(s-a)-2/3r^2 = 0
% INPUT:
%       s    -> Spannungsdeviator
%       a    -> Backstresstensor
%       r    -> Radius FF
%       ds   -> Deviatorisches Spannungsinkrement
%   n_normal -> Anzahl normalenkomponenten
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
if ntens < n_normal
    msg = 'Tensor zu klein oder angegebne Anzahl der Normalenkomponenten zu groß';
    error(msg)
elseif size(a,1) ~= ntens || size(ds,1) ~= ntens
    msg = 'Tensoren haben unterschiedliche Dimension';
    error(msg)
end
devtol = 1e-4;
dummys = sum(s(1:n_normal));
dummya = sum(a(1:n_normal));
if ntens == 6 % 3d Spannungszustand
    if abs(dummys) > devtol 
        msg = 'Übergebenen Spannung ist kein Deviator';
        error(msg)
    elseif abs(dummya) > devtol
        msg = 'Übergebener backstress ist kein Deviator';
        error(msg)
    end
end

% sums = sum(ds);
% if sums == 0
%     xel = 1;
%     return
% end
% Abfangen von Nullen
sums = sum(abs(ds) <= 1e-13);
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
if nargin > 5
    if strcmp(varargin{1,1},'jiang')
        r = sqrt(3) * r;
    elseif strcmp(varargin{1,1},'doring')
        r = sqrt(3/2) * r;
    end
end

% Auswerten 
if ntens == 6 % 3d
    % Notation:
    dummy_rel = [s - a;s(n_normal+1:end)-a(n_normal+1:end)];
    dummy_ds = [ds;ds(n_normal+1:end)];


    q1 = 2 * dummy_rel' * dummy_ds /(dummy_ds' * dummy_ds);
    q2 = (dummy_rel' * dummy_rel - 2/3 * r^2 )/(dummy_ds' * dummy_ds);


elseif ntens == 3 % ESZ
    dummy_rel = [s(1)-a(1);...
                 s(2)-a(2);...
                 -s(1)-s(2)+a(1)+a(2);...
                 s(3)-a(3);...
                 0;...
                 0;...
                 s(3)-a(3);...
                 0;...
                 0];
    dummy_ds = [ds(1);...
                ds(2);...
                -ds(1)-ds(2);...
                ds(3);...
                0;...
                0;...
                ds(3);...
                0;...
                0];
            
    q1 = 2 * dummy_rel' * dummy_ds /(dummy_ds' * dummy_ds);
    q2 = (dummy_rel' * dummy_rel - 2/3 * r^2 )/(dummy_ds' * dummy_ds);

    % fprintf('xel=%.5d\n',xel)
elseif ntens == 2 % sig-tau
    dummy_rel = [s(1)-a(1);...
                 -0.5*(s(1)-a(1));...
                 -0.5*(s(1)-a(1));...
                 s(2)-a(2);...
                 0;...
                 0;...
                 s(2)-a(2);...
                 0;...
                 0];
    dummy_ds = [ds(1);...
                -0.5*ds(1);...
                -0.5*ds(1);...
                ds(2);...
                0;...
                0;...
                ds(2);...
                0;...
                0];
            
    q1 = 2 * dummy_rel' * dummy_ds /(dummy_ds' * dummy_ds);
    q2 = (dummy_rel' * dummy_rel - 2/3 * r^2 )/(dummy_ds' * dummy_ds);

    % fprintf('xel=%.5d\n',xel)
else % einachsige vergleichszustände
    q1 = 2 * (s-a) * ds /ds^2;
    q2 = ((s-a)^2 - r^2 )/(ds^2);

    
end
% Abfangen kleiner Überspannungen
if q2 > 0 && q2 < 1e-6
    q2 = 0;
end
% Berechne elastioschen Anteil
x1 = -q1/2 + sqrt(q1^2/4 - q2);
x2 = -q1/2 - sqrt(q1^2/4 - q2);
xel = max(x1,x2);

% Abbruch bei Fehler
xtol = 1e-4; % setze Toleranz
if xel < -xtol 
    % Mögliche Ursache ist das Konsistenzbedingung nicht eingehalten wurde
    msg='elastischer Anteil kleiner Null berechnet';
    error(msg)
elseif abs(imag(xel)) > 0
    msg=['imaginärer elastischer Teilschritt, ',...
         'Ursache ist zu große Überspannung ', 'q1 = ', num2str(q1), ...
         ';  q2 = ' , num2str(q2), '. q2 ist immer kleiner gleich Null',...
         ' solange keine Überspannung auftritt'];
    error(msg)
elseif xel > 1 + xtol
    msg='elastischer Anteil größer 1 berechnet';
    error(msg)
end
 