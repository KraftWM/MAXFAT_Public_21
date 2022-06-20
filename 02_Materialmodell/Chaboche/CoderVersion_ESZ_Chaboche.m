%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Matlab Coder Version von:                                            %
%    Modellgleichung ESZ Chaboche Modell                                  %
%                                                                         %
%    Aufgerufen in:                                                       %
%    chaboche.m                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dX = ...
          CoderVersion_ESZ_Chaboche(~, X, M, ink, ink_flag, parameter, C, D, ...
                            P, P_hat, A, P_check)
% Konkretes Materialmodell für ESZ Spannungszustände, gibt bei vorgabe eines
% Lastinkrementes die Inkremente der inneren Variablen zurück. Dabei wird
% angenommen, das jedes übergebene Inkrement elastisch-plastische 
% Deformationen hervorruft.
%
% Materialmodell:
% kinematische Verfestigung -> durch Summe von Teilbackstresstensoren
% jeder TBST verhält sich nach Armstrong Frederick Kinematik
%   INPUT:
%         t         -> Zeit (auch wenn nicht explizit gebraucht)
%         X         -> Zustand
%         M         -> Anzahl der Backstresstensoren (auf max 20 begrenzt)
%         ink       -> Lastinkrement
%         ink_flag  -> Belastungsart
%         parameter -> Modellparameter
%                      = 2 + num_alpha * 2 + 1
%         C,D       -> Elastische Steifigkeit und Nachgiebigkeit   
%   P, P_line ...   -> Diverse Abbildungen
%
%
%   OUTPUT:
%         dX -> Inkremente der Zustandsvariablen
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik             |
%  ----------------------------------------------------------------------


%--------------------------------------------------------------------------
%                   Materialparameter
%--------------------------------------------------------------------------

ntens = 3;                                                                 % Tensorkomponenten
% M = (length(parameter)-5)/2;                                               % anzahl Backstresstensoren

% isotrope Verfestigung
q = parameter(3);
gamma = parameter(4);
% kinematische Verfestigung
zeta_i = parameter(5:4+M);
r_i = parameter(5+M:end-1);
r0 = parameter(end);                                                       % startradius fliessfläche
h_i = zeta_i.*r_i;

% oft verwendete Konstanden
w3d2 = sqrt(3/2);
w2d3 = sqrt(2/3);

%--------------------------------------------------------------------------
%                   Zustandsvariablen
%--------------------------------------------------------------------------

% spannungen und dehnungen
if ink_flag == 0 % Spansteu.
    eps = X(1:ntens);
    epsp = X(ntens+1:2*ntens);
    sig = C * (eps - epsp);
elseif ink_flag == 1 % Dehnsteu
    sig = X(1:ntens);
else
    error('Faslche Steuerung angegeben. (0) Spansteu. oder (1) Dehnsteu.')
end
% backstress
alpha = reshape(X(2*ntens+1:(M+2)*ntens),ntens,M);
% radius fließfläche
r = X(end-1);

%--------------------------------------------------------------------------
%                   Normale an die Fließfläche
%--------------------------------------------------------------------------

% Spannungsdeviator
s = P * sig;
% Backstress
a = sum(alpha,2);
% normale an Fließfläche
n = (w3d2/r).* P_hat * (s-a);
nTrans = n.';

%--------------------------------------------------------------------------
%                   Ableitung der Teilbackstresstensoren
%--------------------------------------------------------------------------

% Ableitung Teilbackstresstensoren
dalpha_dp = zeros(ntens,M);
var = w2d3 .* A * n;
for ii = 1 : M
    dalpha_dp(:,ii) = h_i(ii) * var - zeta_i(ii) * alpha(:,ii);
end
% dalpha_dp = w2d3 .* A * n * h_i - zeta_i.*alpha;
% Ableitung gesamtbackstresstensor
da_dp = sum(dalpha_dp,2);

%--------------------------------------------------------------------------
%                   Ableitung des radius der Fließfläche
%--------------------------------------------------------------------------

dr_dp = (q-gamma*(r-r0));

%--------------------------------------------------------------------------
%                   plastischer Tangentenmodul
%--------------------------------------------------------------------------

h = w3d2 * nTrans * P_check * da_dp + dr_dp;

%--------------------------------------------------------------------------
%                   inkremente plastische bogenlänge
%--------------------------------------------------------------------------

if ink_flag == 0 % spannungssteuerung
    dp = (w3d2/h) * (nTrans * ink); 
elseif ink_flag == 1 % dehnungssteuerung
    dp = (w3d2 .* nTrans * (C * ink)) / ...
         (h + 3/2.* nTrans * ( C * n));
else
    error('Faslche Steuerung angegeben. (0) Spansteu. oder (1) Dehnsteu.')
end

%--------------------------------------------------------------------------
%                   inkremente der zustandsvariablen
%--------------------------------------------------------------------------

depsp=dp.*w3d2.*n;
dalpha=dalpha_dp.*dp;
dr=dr_dp * dp;

%--------------------------------------------------------------------------
%                   Zusammenfassen der Inkremente
%--------------------------------------------------------------------------

if ink_flag == 0 % spansteu
    dout = D * ink + depsp;
elseif ink_flag == 1 % dehnsteu
    dout = C * (ink - depsp);
else
    error('Faslche Steuerung angegeben. (0) Spansteu. oder (1) Dehnsteu.')
end
dX = [dout; depsp ; reshape(dalpha,ntens*M,1); dr; dp]; 
    
end % Ende Modell ESZ