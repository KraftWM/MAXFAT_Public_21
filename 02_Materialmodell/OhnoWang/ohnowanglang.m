function [ZVARneu] = ohnowanglang(desig,ZVAR,para,epara)
% Ohno Wang Plastizitätsmodell für pseudo Stress approach nach Asatz 
% von Lang, implementierung nach Döring 
%
% ! Nur Spannungssteuerung im ESZ
%
% QUELLE:
%       Ohno et al. 1993 -  KINEMATIC HARDENING RULES WITH CRITICAL
%                           STATE OF DYNAMIC RECOVERY, PART I
%       Lang et al       - A multiaxial stress-strain correction scheme
%
% INPUT:
%  dESIG -> Inkrement in pseudo Spannungen
%  ZVAR  -> Zustandsvariablen
%  para  -> Parameter des Materialmodells
% epara  -> Parameter des Strukturmodells
%
% OUTPUT:
%  ZVARneu -> neuer Zustand nach Inkrement
%
%__________________________________________________________________________
%
% Zustandsvariablen:
% Zuerst Zustandsvariablen des Struckturmodells (wie bei Spannungssteuerung)
% dann zusätzliche Variablen des Materialmodells 
% ZVAR = [eeps    -> Pseudo Dehnungen (! nicht pseudo elastische Dehnungen)
%         epsp    -> "reale" Plastische Dehnungen 
%         ealphai -> pseudo Backstress
%         p       -> "reale" plastische Bogenlänge
%         alphai  -> "reale" Backstresstensoren
% eeps = ZVAR(1:ntens)
% epsp = ZVAR(ntens+1:2*ntens)
% ealphai = ZVAR(2*ntens+1:(2+eM)*ntens)
% p = ZVAR((2+eM)*ntens+1)
% alphai = ZVAR((2+eM)*ntens+2:(2+eM+M)*ntens+1)
%
% Darstellung von Tensoren
%         sig11              eps11
%  sig =  sig22       eps =  eps22
%         sig12             2eps12
%__________________________________________________________________________          
%
% Parameter:
% zuerst Material- dann Struckturmodell 
%     para = [E, nu,  c_i,  r_i,  chi_i,  r0]
%       M = (length(para)-3)/3;
%     epara = [E, nu,  ec_i, er_i, echi_i, er0]
%       eM = (length(epara)-3)/3; 
%
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft                                                   |
% |  Stand: Januar 2020                                                 |
%  ----------------------------------------------------------------------

% -------------------------------------------------------------------------
% !!!!! AKTUELL NUR ESZ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Für weitere Spannungszustände müssten "nur" die Modellgleichung
% implementiert werden
% -------------------------------------------------------------------------
ntens = 3;
ndi = 2;


%--------------------------------------------------------------------------
%                  Materialparamter ermitteln                                  
%--------------------------------------------------------------------------
M = (length(para)-3)/3;                                                    % Anzahl Backstresstensoren
eM = (length(epara)-3)/3;                                                  % Anzahl Backstresstensoren Strukturmodell

E = para(1);                                                               % Elastizitätsmodul
nu = para(2);                                                              % Querdehnzahl   
er0 = epara(3*eM+3);                                                       % Radius pseudo Fließfläche


%--------------------------------------------------------------------------
%                        Hilfsmatrizen                                    
%--------------------------------------------------------------------------

% statische Matrizen
[P, P_line] = set_maps(ntens,ndi);
% Steifigkeit
CEL = elast_steifigkeit(E,nu,ntens,ndi);
% Nachgiebigkeit
DEL = elast_nachgiebigkeit(E,nu,ntens,ndi);

%--------------------------------------------------------------------------
%            Identifikation der Zustandsvariablen                         
%--------------------------------------------------------------------------

% pseudo gesamtdehnung 
eeps =  ZVAR(1:ntens);
% plastische dehnung
epsp = ZVAR(ntens+1:2*ntens);
% Pseudo (elastische9 Spannungen
esig = CEL * (eeps - epsp);
% pseudo backstress
ealphai = reshape(ZVAR(2*ntens+1:(eM+2)*ntens),ntens,eM);

%--------------------------------------------------------------------------
%               elastischer Trial Step                         
%--------------------------------------------------------------------------

es = P * esig;                                                             % Pseudo Spannungsdeviator
ea = sum(ealphai,2);                                                       % Pseudo Gesamtbackstress
des = P * desig;                                                           % Inkrement Pseudo Spannungsdeviator
es_tr = es + des;                                                          % Pseudo Versuchsspannungsdeviator
ebeta = es_tr - ea;                                                        % Pseudo effektive Spannung
if ntens == 1 % 1D
    F_tr = abs(ebeta) - er0;                                               % Pseudo Überspannung
else % Alles andere
    F_tr = ebeta' * P_line * ebeta - 2/3 * er0^2;                          % Pseudo Überspannung
end
FTOL = 1e-7;

%--------------------------------------------------------------------------
%                   Trial Step wird angenommen
%--------------------------------------------------------------------------
if F_tr < FTOL

    % Updaten der pseudo gesamtdehnung
    eeps = eeps + DEL * desig;      
    % Updaten Zustandsvariablen
    ZVARneu = [eeps;ZVAR(ntens+1:end)];

%--------------------------------------------------------------------------
%                   Trial Step wird nicht angenommen
%--------------------------------------------------------------------------
else
    
    % elastischer trial step
%     xel = elastink(es,ea,er0,des,ndi);                                     % Elastischer Anteil des Inkrements
    xel = elastink2(es,ea,P_line,er0,des,FTOL);
    % elastischen Anteil aufbringen
    eeps = eeps + DEL * (xel*desig);                                       % Pseudo Dehnungen
    ZVAR(1:ntens) = eeps;
    
    % restinkrement
    desig = (1-xel) * desig;

    % Integration des plastischen Anteils
    options = [];
    if ntens == 6 % 3D
        
        msg = 'nicht implementiert';
        error(msg)
        
    elseif ntens == 3 && ndi == 2 % ESZ
        
        % Abbildungen
        [P, P_line, P_hat, A, P_check] = set_maps(ntens,ndi);
        % Integration
        
        % Matlab Version
%         [~,ZVARneu] = rk87(@modellESZ,[0,1], ZVAR,options,...
%                             desig, para, epara, CEL, DEL, ...
%                             P, P_line, P_hat, A, P_check);

        % Mes Funktion
        [~,ZVARneu] = rk87(@CoderVersion_ESZ_OhnoWangLang_mex,[0,1], ZVAR,options,...
                            M, eM, desig, para, epara, CEL, DEL, ...
                            P, P_line, P_hat, A, P_check);
                        
    elseif ntens == 1 % 1D
        
        msg = 'nicht implementiert';
        error(msg)
        
    end
    
    % nur letzten Schritt ausgeben
    ZVARneu = ZVARneu(end,:)'; 
    
    % Prüfe die Konsistenzbedingung
    ZVARneu = konsistenzbedingung(ZVARneu,eM,P,P_line,ntens,er0,CEL,DEL);
       
end
end % ENde Hauptfunktion



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Prüfe die konsistenzbedingung                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = konsistenzbedingung(X,eM,P,P_line,ntens,r,C,D)
% Window Ausgabe
% fprintf('Zustandsvariablen Original:\n')
% fprintf('%.32d \n',X)

% auslesen zustände
% Pseudo Dehnungen (! nicht pseudo elastische Dehnungen)
eeps = X(1:ntens);
% "reale" Plastische Dehnungen 
epsp = X(ntens+1:2*ntens);
% pseudo spannungen 
esig = C * (eeps - epsp);
% pseudo Backstress
ealphai = reshape(X(2*ntens+1 : (eM+2)*ntens),3,eM);
ea = sum(ealphai,2);
% "reale" Backstresstensoren
% alphai = reshape(X((M+2)*ntens+2 : (2*M+2)*ntens+1),3,M);


% Relativspannung
es = P * esig;
ebeta = es - ea;
% fprintf('Spanndev orig:\n')
% fprintf('%.32d\n',s)
% Fließfläche
if ntens == 1
    SV = ebeta^2;
else
    SV = 1.5 * ebeta' * P_line * ebeta;
end
% F1 = SV - r^2;
% fprintf('Überspannung Orig: %.4d \n',F1)

% Prüfe Konsistenzbedingung
if SV - r^2 > 1e-11
    
    % Fehler/Warnung wenn zu stark verletzt
    if SV - r^2 > 1e-1
        msg = 'Konsistenzbedinung verletzt';
        error(msg)
    elseif SV - r^2 > 1e-4
        msg = 'Konsistenzbedinung geringfügig verletzt';
        warning(msg)
    end

    % Korrigiere Spannungsdeviator (Rel. Span. radial auf FF zurückprojezieren)
    fak = r/sqrt(SV);
    
    % ... ESZ nur Deviator korriegieren
    es = ea + fak*ebeta;
    % korrigiere Spannungen
    esig = [2 1 0; 1 2 0; 0 0 1] * es;


    % Korrigierte Zustandsgrößen
    eepse = D*esig;
    eeps = eepse+epsp;
    X(1:ntens) = eeps;

end
end % Ende Prüfen Konsitenzbedingung





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Modellgleichung ESZ                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dX = modellESZ(~, X, desig, para, epara, CEL, DEL, ...
                        P, P_line, P_hat, A, P_check)
% Inkrements des System von DGL fürs Ohno Wang Modell mit den Ansatz von Lang
%
% INPUT:
%      ~       -> Zeit wird hier nicht gebraucht, aber wird von ode solver
%                 benötigt
%     X        -> aktueller Zustand
%  para        -> Material- und Struckturparameter
%  desig       -> inkrement des pseudospannungstensors
% CEL,DEL      -> Elastische Steifigkeit und Nachgiebigkeit   
% P, P_line... -> Diverse Abbildungen
%
% OUTPUT:
%   dX   -> Inkrement der Zustandsgrößen
%
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
%    Parameter
% -------------------------------------------------------------------------

% Spannungszustand
ntens = 3;
% ndi = 2;
% Anzahl Backstresstensoren
M = (length(para)-3)/3;
eM = (length(epara)-3)/3;                                                    % Anzahl Backstresstensoren
   
% Material
c_i = para(3:2+M);
r_i = para(3+M:2+2*M);
chi_i = para(3+2*M:2+3*M);
% Strucktur
ec_i = epara(3:2+eM);
er_i = epara(3+eM:2+2*eM);
echi_i = epara(3+2*eM:2+3*eM);
er0 = epara(3+3*eM);


%--------------------------------------------------------------------------
%            Zustandsvariablen
%--------------------------------------------------------------------------
% Pseudo Dehnungen (! nicht pseudo elastische Dehnungen)
eeps = X(1:ntens);
% "reale" Plastische Dehnungen 
epsp = X(ntens+1:2*ntens);
% pseudo spannungen 
esig = CEL * (eeps - epsp);
% pseudo Backstress
ealphai = reshape(X(2*ntens+1 : (eM+2)*ntens),3,eM);
% "reale" Backstresstensoren
alphai = reshape(X((eM+2)*ntens+2 : (eM+M+2)*ntens+1),3,M);
% numerisch null
delta = 1e-40;

%--------------------------------------------------------------------------
%               Normale an die Strucktur Fließfläche
%--------------------------------------------------------------------------

% pseudobackstress
ea = sum(ealphai,2);
% aktueller pseudo Spannungsdeviator
es = P * esig;
% pseudo effektive spannung
ebeta = es  - ea;
% Normale
en = (sqrt(3/2)/er0) * P_hat * ebeta;
transen = en';

%--------------------------------------------------------------------------
%               Ableitung der Teilbackstresstensoren
%--------------------------------------------------------------------------

% normen der teilbackstresstensoren
norm_ai = sqrt( sum( (P_line * alphai) .* alphai ) );
norm_ai(norm_ai == 0) = delta;

% normen der pseudo teilbackstresstensoren
norm_eai = sqrt( sum( (P_line * ealphai) .* ealphai ) );
norm_eai(norm_eai == 0) = delta;

% hilfsgrößen
Li = alphai./norm_ai;
eLi = ealphai./norm_eai;

% % weitere Hilfsgrößen
% dummy1 = norm_ai./r_i;                                                     % Nähe Backstress zu Begrenzungsfläche
% dummy1( dummy1 > 1 ) = 1;
% dummy1 = dummy1.^chi_i;
% edummy1 = norm_eai./er_i;                                                  
% edummy1( edummy1 > 1 ) = 1;
% edummy1 = edummy1.^echi_i;
% 
% dummy2 = transen * P_check * Li;                                          
% dummy2 = 0.5 * (dummy2 + abs(dummy2));
% edummy2 = transen * P_check * eLi; 
% edummy2 = 0.5 * (edummy2 + abs(edummy2));
% edummy2 = edummy1 .* edummy2 .* ealphai;
% dummy3 = A * en * r_i;
% edummy3 = A * en * er_i;
% 
% % Ableitungen teilbackstress
% dalphai_dp = c_i .* (dummy3 - dummy1 .* dummy2 .* alphai); 
% 
% % Ableitung pseudo teilbackstress
% dealphai_dp = ec_i .* (edummy3 - edummy2); 

% weitere hilfsgrößen ( !!!! IMPLEMENTIERUNG DER SCHLEIFE SCHEINT SCHNELLER
% ZU SEIN !!!! )
var1 = c_i .* r_i;
evar1 = ec_i .* er_i;
var2 = A * en;
var3 = norm_ai./r_i;
var3(var3>1) = 1;
evar3 = norm_eai./er_i;
evar3(evar3>1) = 1;
var4 = transen * P_check * Li;
var4 = 0.5 * (var4 + abs(var4));
evar4 = transen * P_check * eLi;
evar4 = 0.5 * (evar4 + abs(evar4));

% Schleife über alle Backstresstensoren
% init ableitungen
dealphai_dp = zeros(ntens,eM);
dalphai_dp = zeros(ntens,M);
for ii = 1 : M
    
    % backstress
    dalphai_dp(:,ii) = var1(ii) * (var2 - (var3(ii)).^(chi_i(ii)+1).*var4(ii).*Li(:,ii));
        
end

for ii = 1 : eM
       
    % pseudo backstress
    dealphai_dp(:,ii) = evar1(ii) * ( var2 - (evar3(ii)).^(echi_i(ii)+1).*evar4(ii).*eLi(:,ii) );
    
end

% Ableitung gesamter pseudo backstress
dea_dp = sum(dealphai_dp,2);

%--------------------------------------------------------------------------
%                  plastischer tangentenmodul                   
%--------------------------------------------------------------------------

eh = transen * P_check * dea_dp;

%--------------------------------------------------------------------------
%                  inkrement plastische Bogenlänge                   
%--------------------------------------------------------------------------

dp = (transen * desig)/eh;

%--------------------------------------------------------------------------
%                  inkremente der zustandsvariablen                  
%--------------------------------------------------------------------------
depsp = en .* dp;
deeps = depsp + DEL * desig;
dalphai = dalphai_dp .* dp;
dealphai = dealphai_dp .* dp;

%--------------------------------------------------------------------------
%                   Zusammenfassen der Inkremente                        
%-------------------------------------------------------------------------- 
  
dX = [deeps; depsp; reshape(dealphai,3*eM,1); dp; reshape(dalphai,3*M,1)];
end
