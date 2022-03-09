%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Matlab Coder Version von:                                            %
%    Modellgleichung ESZ Chaboche Modell f¸r Integration mit              %
%    Energieinkrement (nach mir)                                          %
%                                                                         %
%    Aufgerufen in:                                                       %
%    chabocheenergyexp.m                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dX = CoderVersion_ESZ_ChabocheEnergie(~, ZVAR, domega, para, ...
                       GINV, PHAT,...
                       M_, MHAT, A, MCHECK)
% Inkrement des Chaboche Modells bei energiegesteuerter Integration                   
%
% INPUT:
%      ~       -> Zeit wird hier nicht gebraucht, aber wird von ode solver
%                 benˆtigt
%  ZVAR        -> aktueller Zustand
%  para        -> Material- und Struckturparameter
%  domega      -> inkrement der energien
%   GINV       -> inverse Ableitung Energie nach Spannungen
%  PHAT        -> ABleitung Energie nach plastischen Dehnungen (!als vektor)
% M_, MLINE... -> Diverse Abbildungen
%
% OUTPUT:
%   dX   -> Inkrement der Zustandsgrˆﬂen                   
%
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
%                   Materialparameter
%--------------------------------------------------------------------------

% Zuweisen der Materialparameter
ntens = 3;                                                                 % Tensorkomponenten
M = (length(para)-5)/2;                                               % anzahl Backstresstensoren

r0 = para(end);                                                       % startradius fliessfl‰che
q = 0;
gamma = 0;

% kinematische Verfestigung
zeta_i = para(5:4+M);
r_i = para(5+M:end-1);
r0 = para(end);                                                       % startradius fliessfl‰che
h_i = zeta_i.*r_i;

% oft verwendete Konstanden
w3d2 = sqrt(3/2);
w2d3 = sqrt(2/3);

%--------------------------------------------------------------------------
%                   Zustandsvariablen
%--------------------------------------------------------------------------

% Spannungen
sig = ZVAR(1:ntens);
% backstresstensoren
alpha = reshape(ZVAR(2*ntens+1:(M+2)*ntens),ntens,M);
% radius flieﬂfl‰che
r = ZVAR(end-1);

%--------------------------------------------------------------------------
%                   Normale an die Flieﬂfl‰che
%--------------------------------------------------------------------------

% Spannungsdeviator
s = M_ * sig;
% Backstress
a = sum(alpha,2);
% normale an Flieﬂfl‰che
n = (w3d2/r0).* MHAT * (s-a);
transn = n.';

%--------------------------------------------------------------------------
%                   Ableitung der Teilbackstresstensoren
%--------------------------------------------------------------------------

% Ableitung Teilbackstresstensoren
dalpha_dp = zeros(ntens,M);
dummy = w2d3 * A * n;
for ii = 1:M
    dalpha_dp(:,ii) = dummy * h_i(ii) - zeta_i(ii) * alpha(:,ii);
end
% Ableitung gesamtbackstresstensor
da_dp = sum(dalpha_dp,2);

%--------------------------------------------------------------------------
%                   Ableitung des radius der Flieﬂfl‰che
%--------------------------------------------------------------------------

dr_dp = (q-gamma*(r-r0));


%--------------------------------------------------------------------------
%                   plastischer Tangentenmodul
%--------------------------------------------------------------------------

h = w3d2 * transn * MCHECK * da_dp + dr_dp;

%--------------------------------------------------------------------------
%                   inkremente plastische bogenl‰nge
%--------------------------------------------------------------------------
dummy1 = GINV *domega;
if size(PHAT,2) == 1
    dummy2 = GINV * ( PHAT .* n );
else
    dummy2 = GINV * ( PHAT * n );
end
dp = (w3d2*transn * dummy1)/(h + 3/2 * transn * dummy2 );

%--------------------------------------------------------------------------
%                   inkremente der zustandsvariablen
%--------------------------------------------------------------------------

depsp=dp.*w3d2.*n;
dalpha=dalpha_dp.*dp;
dr=dr_dp * dp;
dsig = dummy1 - dummy2.*(w3d2*dp);


%--------------------------------------------------------------------------
%                   Zusammenfassen der Inkremente
%--------------------------------------------------------------------------

% dX = [dsig; depsp ; reshape(dalpha,ntens*M,1); dp]; 
dX = [dsig; depsp ; dalpha(:); dr; dp]; 

end % Ende Modellgleichung ESZ