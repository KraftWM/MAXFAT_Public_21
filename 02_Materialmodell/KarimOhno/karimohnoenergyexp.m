function [ZVARneu] = karimohnoenergyexp(domega,ZVAR,para,GHAT,PHAT)
% Karim Ohno Plastizit�tsmodell.
% implementiert f�r die energiegesteuerte integration bei inkrementellen
% Kerbn�herungen
%
% INPUT:
%  domega -> Inkrement der energie
%  ZVAR   -> Zustandsvariablen
%  para   -> Parameter des Pseudo Modells und des Materialmodells
%  CEL    -> elastische Steifigkeit
%  GHAT   -> Ableitung energie nach spannungen
%  PHAT   -> Ableitung energie nach plastischen Dehnungen
%
% OUTPUT:
%  ZVARneu -> neuer zustand nach inkrement
%
%__________________________________________________________________________
%
% Zustandsvariablen:
% Zuerst so wie bei dehnungsgesteuerter integration 
% ZVAR = [sig;epsp;alphai;p]
%
% Darstellung von Tensoren
%         sig11              eps11
%  sig =  sig22       eps =  eps22
%         sig12             2eps12
%__________________________________________________________________________          
%
% Parameter:
% Materialparameter des Modells
%                      [E,nu,c_i,r_i,mu_i,r0]
%                       mu_i(i) = 0 -> OhnoWang Modell 1 Kinematik
%                       mu_i(i) = 1 -> Armstrong Frederick Kinematik
%
%       M = (length(para) - 3)/3
%
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft                                                   |
% |  Stand: Juli 2020                                                 |
%  ----------------------------------------------------------------------

% -------------------------------------------------------------------------
% !!!!! AKTUELL NUR ESZ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% F�r weitere Spannungszust�nde m�ssten "nur" die Modellgleichung
% implementiert werden
% -------------------------------------------------------------------------
ntens = 3;
ndi = 2;

%--------------------------------------------------------------------------
%                  Materialparamter ermitteln                                  
%--------------------------------------------------------------------------

% Anzahl Backstress
M = (length(para)-3)/3;
% Elastizit�tskonstanten
r0 = para(end);

%--------------------------------------------------------------------------
%                        Hilfsmatrizen                                    
%--------------------------------------------------------------------------
% Abbildungen 
[M_, MLINE] = set_maps(ntens,ndi);                                                                           
% inverse ableitung nach Spannungen
dummy = GHAT(1,1)*GHAT(2,2)-GHAT(1,2)*GHAT(2,1);
GINV = (1/dummy) .* [ GHAT(2,2),-GHAT(1,2),0;...
                     -GHAT(2,1), GHAT(1,1),0;...
                              0,         0, dummy/GHAT(3,3)];
                          
%--------------------------------------------------------------------------
%                  Zustandsvariablen auslesen                                 
%--------------------------------------------------------------------------

% Spannungen und Dehnungen
sig = ZVAR(1:ntens);
% Backstresstensoren
alpha = reshape( ZVAR(2*ntens+1:(M+2)*ntens) , ntens, M);

%--------------------------------------------------------------------------
%                       Elastischer Trial Step
%--------------------------------------------------------------------------
dsig_tr = GINV * domega;                                                   % elastisches Versuchsinkrement
s = M_ * sig;                                                              % Span dev
ds = M_ * dsig_tr;                                                      % Versuchsinkrement Spandeviator 
s_tr = s + ds;
a = sum(alpha,2);                                                          % Gesamtbackstress
beta = s_tr - a;                                                           % trial effektive spannung
F_tr = beta' * MLINE * beta - 2/3 * r0^2;                                 % Trial Flie�funktion
FTOL = 1e-7;                                                              % Toleranz f�r Abweichungen F ~= 0

%--------------------------------------------------------------------------
%                   Trial Step wird angenommen
%--------------------------------------------------------------------------
if F_tr < FTOL
    
    % elastisches Update
    sig = sig + dsig_tr;
    ZVARneu = [sig ; ZVAR(ntens+1:end)];

%--------------------------------------------------------------------------
%                   Trial Step wird abgelehnt
%--------------------------------------------------------------------------
else
   
    % 1. elastischer Anteil des Inkrements
%     xel = elastink(s,a,r0,ds,ndi);
    xel = elastink2(s,a,MLINE,r0,ds,FTOL);  
    % 2. elastisches update
    sig = sig + xel * dsig_tr;
    ZVAR = [sig ; ZVAR(ntens+1:end)];
    
    % 3. restinkrement
    domega = (1-xel)*domega;
    
    % 4. Integration des restinkrements
    options = [];
    if ntens == 6 % 3D
        
        msg = 'nicht implementiert';
        error(msg)
        
    elseif ntens == 3 && ndi == 2 % ESZ
        
        % Abbildungen
        [M_, MLINE, MHAT, A, MCHECK] = set_maps(ntens,ndi);
        % Integration
        
        % Matlab Version
        [~,ZVARneu] = rk87(@modellESZ,[0,1], ZVAR,options,...
                            domega, para, GINV, PHAT, ...
                            M_, MLINE, MHAT, A, MCHECK);

                        
    elseif ntens == 1 % 1D
        
        msg = 'nicht implementiert';
        error(msg)
        
    end
    
    % nur letzten Schritt ausgeben
    ZVARneu = ZVARneu(end,:)'; 
    
end % Ende Verzweigung elastischer Trial Step
end % Ende Funktion 














%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Modellgleichungen ESZ                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dX = modellESZ(~, ZVAR, domega, para, ...
                       GINV, PHAT,...
                       M_, MLINE, MHAT, A, MCHECK)
% Inkrement des Karim Ohno Modells bei energiegesteuerter Integration                   
%
% INPUT:
%      ~       -> Zeit wird hier nicht gebraucht, aber wird von ode solver
%                 ben�tigt
%  ZVAR        -> aktueller Zustand
%  para        -> Material- und Struckturparameter
%  domega      -> inkrement der energien
%   GINV       -> inverse Ableitung Energie nach Spannungen
%  PHAT        -> ABleitung Energie nach plastischen Dehnungen (!als vektor)
% M_, MLINE... -> Diverse Abbildungen
%
% OUTPUT:
%   dX   -> Inkrement der Zustandsgr��en                   
%
% -------------------------------------------------------------------------

%--------------------------------------------------------------------------
%                   Materialparameter
%--------------------------------------------------------------------------

% Zuweisen der Materialparameter
ntens = 3;                                                                 % Tensorkomponenten
M = (length(para)-3)/3;                                               % anzahl Backstresstensoren

r0 = para(end);                                                       % startradius fliessfl�che

% kinematische Verfestigung
c_i = para(3:2+M);
r_i = para(3+M:2+2*M);
mu_i = para(3+2*M:2+3*M);
h_i = c_i .* r_i;
% oft verwendete Konstanden
w2d3 = sqrt(2/3);
w3d2 = sqrt(3/2);

delta = 1e-40;                                                             % numerisch 0

%--------------------------------------------------------------------------
%                   Zustandsvariablen
%--------------------------------------------------------------------------

% Spannungen
sig = ZVAR(1:ntens);
% backstresstensoren
alpha = reshape(ZVAR(2*ntens+1:(M+2)*ntens),ntens,M);

%--------------------------------------------------------------------------
%                   Normale an die Flie�fl�che
%--------------------------------------------------------------------------

% Spannungsdeviator
s = M_ * sig;
% Backstress
a = sum(alpha,2);
% normale an Flie�fl�che
n = (w3d2/r0).* MHAT * (s-a);
transn = n.';

%--------------------------------------------------------------------------
%                   Ableitung der Teilbackstresstensoren
%--------------------------------------------------------------------------

% Normen der Teilbackstresstensoren
norm_ai = sqrt( sum( (MLINE*alpha) .* alpha ) );
norm_ai(norm_ai == 0) = delta; 

% Hilfvariablen
% Li = alpha./norm_ai;
% var1 = c_i .* r_i;
Li = w2d3*alpha./norm_ai;
var2 = A * n;
var3 = w3d2*norm_ai./r_i;
var3(var3 >= 1) = 1;
var3(var3 < 1) = 0;
var4 = w3d2 * transn * MCHECK * Li - mu_i;
var4 = 0.5 * (var4 + abs(var4));

% Ableitung Teilbackstresstensoren
dalpha_dp = zeros(ntens,M);
for ii = 1 : M

    dalpha_dp(:,ii) = h_i(ii) * ( ...
                                w2d3 * var2 - ...
                                mu_i(ii)/r_i(ii) .* alpha(:,ii) - ...
                                var3(ii) * var4(ii) .* Li(:,ii) ...
                                );
    
end
% Ableitung gesamtbackstresstensor
da_dp = sum(dalpha_dp,2);

%--------------------------------------------------------------------------
%                   plastischer Tangentenmodul
%--------------------------------------------------------------------------

h = w3d2 * transn * MCHECK * da_dp;

%--------------------------------------------------------------------------
%                   inkremente plastische bogenl�nge
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
dsig = dummy1 - dummy2.*(w3d2*dp);

%--------------------------------------------------------------------------
%                   Zusammenfassen der Inkremente
%--------------------------------------------------------------------------

% dX = [dsig; depsp ; reshape(dalpha,ntens*M,1); dp]; 
dX = [dsig; depsp ; dalpha(:); dp]; 

end % Ende Modellgleichung ESZ