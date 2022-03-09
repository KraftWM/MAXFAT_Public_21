function [ZVARneu] = jianglang(desig,ZVAR,para,epara)
% Jiang Plastizitätsmodell für pseudo Stress approach nach Asatz von Lang
%
% ! Nur Spannungssteuerung im ESZ
%
% QUELLE:
%       Jiang, Sehitoglu - Modeling of Cyclic Rachtetting Plasticity
%       Lang et al       - A multiaxial stress-strain correction scheme
%
% INPUT:
%  dESIG -> Inkrement in pseudo spannungen
%  ZVAR  -> Zustandsvariablen
%  para  -> Parameter des Pseudo Modells und des Materialmodells
%
% OUTPUT:
%  ZVARneu -> neuer zustand nach inkrement
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
%         erm     -> Pseudo memory surface
%         alphai  -> "reale" Backstresstensoren
%         rm      -> "reale" memory surface
%
% eeps = ZVAR(1:ntens);
% epsp = ZVAR(1+ntens:2*ntens);
% ealphai = ZVAR(2*ntens+1:(2+eM)*ntens);
% p = ZVAR((2+eM)*ntens+1);
% erm = ZVAR((2+eM)*ntens+2);
% alphai = ZVAR((2+eM)*ntens+3:(2+eM+M)*ntens+2);
% rm = ZVAR((2+eM+M)*ntens+3);
%
% Darstellung von Tensoren
%         sig11              eps11
%  sig =  sig22       eps =  eps22
%         sig12             2eps12
%__________________________________________________________________________          
%
% Parameter:
% zuerst Material- dann Struckturmodell 
%     para = [E, nu, a_chi,  b_chi,  ak,  ck,  k1,  cm,  c_i0,  a_i1,  b_i1,
%              a_i2,  b_i2,  r_i,  Q_i]
%   epara = [E, nu, ea_chi, eb_chi, eak, eck, ek1, ecm, ec_i0, ea_i1, eb_i1,
%             ea_i2, eb_i2, er_i, eQ_i]
%       M = (length(para) - 8)/7
%      eM = (length(epara) - 8)/7
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft                                                   |
% |  Stand: Dezember 2019                                               |
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
% M = (length(para)-14)/14;                                                  % Anzahl Backstresstensoren
M = (length(para)-8)/7;                                               % Anzahl Backstresstensoren
eM = (length(epara)-8)/7;                                               % Anzahl Backstresstensoren Strukturmodell
E = para(1);                                                               % Elastizitätsmodul
nu = para(2);                                                              % Querdehnzahl   

% Struckturmodell
eak = epara(5);                                                       % parameter zum bestimmen der Größe der FF
eck = epara(6);
ek1 = epara(7);


%--------------------------------------------------------------------------
%                        Hilfsmatrizen                                    
%--------------------------------------------------------------------------

% statische Matrizen
[P, P_line] = set_maps(3,2);
CEL = elast_steifigkeit(E,nu,ntens,ndi);                                   % Steifigkeit                                       
DEL = elast_nachgiebigkeit(E,nu,ntens,ndi);                                % Nachgiebigkeit

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
% pseudo gedächnissfläche
erm = ZVAR((eM+2)*ntens+2);

%--------------------------------------------------------------------------
%               Radien der Fließflächen berechnen                          
%--------------------------------------------------------------------------

% k = k1 * ( 1 + ak * exp( ck * rm ) );                                      % realer radius
ek = ek1 * ( 1 + eak * exp ( eck * erm ) );                                % pseudo radius

%--------------------------------------------------------------------------
%               elastischer Trial Step                         
%--------------------------------------------------------------------------

es = P * esig;                                                             % Pseudo Spannungsdeviator
ea = sum(ealphai,2);                                                       % Pseudo Gesamtbackstress
des = P * desig;                                                           % Inkrement Pseudo Spannungsdeviator
es_tr = es + des;                                                          % Pseudo Versuchsspannungsdeviator
ebeta = es_tr - ea;                                                        % Pseudo effektive Spannung
if ntens == 1 % 1D
    F_tr = abs(ebeta) - sqrt(3) * ek;                                      % Pseudo Überspannung
else % Alles andere
    F_tr = ebeta' * P_line * ebeta - 2 * ek^2;                             % Pseudo Überspannung
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
%     xel = elastink(es,ea,ek,des,ndi,'jiang');                              % Elastisches Inkrement berechnen
    xel = elastink2(es,ea,P_line,ek,des,FTOL,'jiang');  
    % elastischen Anteil aufbringen
    eeps = eeps + DEL * (xel*desig);
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
        
        % Matlab version
%         [~,ZVARneu] = rk87(@modellESZ,[0,1], ZVAR,options,...
%                             desig, para, epara, CEL, DEL, ...
%                             P, P_line, P_hat, A, P_check);

        % Mex Version
        [~,ZVARneu] = rk87(@CoderVersion_ESZ_JiangLang_mex,[0,1], ZVAR,options,...
                            M, eM, desig, para, epara, CEL, DEL, ...
                            P, P_line, P_hat, A, P_check);
                        
    elseif ntens == 1 % 1D
        
        msg = 'nicht implementiert';
        error(msg)
        
    end
    
    % nur letzten Schritt ausgeben
    ZVARneu = ZVARneu(end,:)'; 
    
    % Testen der Konsistenzbedingung
%     ZVARneu = konsistenzbedingung(ZVARneu,eM,P,P_line,ntens,CEL,DEL,ek1,eak,eck);
       
end
end












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Prüfe die konsistenzbedingung                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = konsistenzbedingung(X,eM,P,P_line,ntens,C,D,ek1,eak,eck)
% Window Ausgabe
% fprintf('Zustandsvariablen Original:\n')
% fprintf('%.32d \n',X)


% Identifikation der Zustandsvariablen
% Pseudo Dehnungen (! nicht pseudo elastische Dehnungen)
eeps = X(1:ntens);
% "reale" Plastische Dehnungen 
epsp = X(ntens+1:2*ntens);
% pseudo spannungen 
esig = C * (eeps - epsp);
% pseudo Backstress
ealphai = reshape(X(2*ntens+1 : (eM+2)*ntens),3,eM);
ea = sum(ealphai,2);  
% Pseudo memory surface
erm = X((eM+2)*ntens+2);

% Abfangen von Nullen 
if erm == 0
    erm = 1e-40;
end
% radius Pseudo Fließfläche
ek = ek1 * ( 1 + eak * exp( eck * erm ) );


% Relativspannung
es = P * esig;
ebeta = es - ea;
% fprintf('Spanndev orig:\n')
% fprintf('%.32d\n',s)
% Fließfläche
if ntens == 1
    SV = 1/3 * ebeta^2;
else
    SV = 0.5 * ebeta' * P_line * ebeta;
end
% F1 = SV - k^2;
% fprintf('Überspannung Orig: %.4d \n',F1)

% Prüfe Konsistenzbedingung
if SV - ek^2 > 1e-11
    
    % Fehler/Warnung wenn zu stark verletzt
    if SV - ek^2 > 1e-1
        msg = 'Konsistenzbedinung verletzt';
        error(msg)
    elseif SV - ek^2 > 1e-4
        msg = 'Konsistenzbedinung geringfügig verletzt';
        warning(msg)
    end

    % Korrigiere Spannungsdeviator (Rel. Span. radial auf FF zurückprojezieren)
    fak = ek/sqrt(SV);
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
% Inkrements des System von DGL fürs Jiang Modell mit den Ansatz von Lang
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
% Anzahl Backstresstensoren
M = (length(para)-8)/7; 
eM = (length(epara)-8)/7;   
% Material
a_chi = para(3);
b_chi = para(4);
cm = para(8);
c_i0 = para(9:9+M-1);
a_i1 = para(9+M:9+2*M-1);
b_i1 = para(9+2*M:9+3*M-1);
a_i2 = para(9+3*M:9+4*M-1);
b_i2 = para(9+4*M:9+5*M-1);
r_i = para(9+5*M:9+6*M-1);
Q_i = para(9+6*M:9+7*M-1);
% Strucktur
ea_chi = epara(3);
eb_chi = epara(4);
eak = epara(5);
eck = epara(6);
ek1 = epara(7);
ecm = epara(8);
ec_i0 = epara(9:9+eM-1);
ea_i1 = epara(9+eM:9+2*eM-1);
eb_i1 = epara(9+2*eM:9+3*eM-1);
ea_i2 = epara(9+3*eM:9+4*eM-1);
eb_i2 = epara(9+4*eM:9+5*eM-1);
er_i = epara(9+5*eM:9+6*eM-1);
eQ_i = epara(9+6*eM:9+7*eM-1);

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
% "reale" plastische Bogenlänge
p = X((eM+2)*ntens+1);
% Pseudo memory surface
erm = X((eM+2)*ntens+2);
% "reale" Backstresstensoren
alphai = reshape(X((eM+2)*ntens+3 : (eM+M+2)*ntens+2),3,M);
% "reale" memory surface
rm = X((eM+M+2)*ntens+3);
% numerisch null
delta = 1e-40;
% Abfangen von Nullen 
if rm == 0
    rm = delta;
end
if erm == 0
    erm = delta;
end
% radius Pseudo Fließfläche
ek = ek1 .* ( 1 + eak * exp( eck * erm ) );

%--------------------------------------------------------------------------
%               Normale an die Strucktur Fließfläche
%--------------------------------------------------------------------------

% pseudobackstress
ea = sum(ealphai,2);
normea = sqrt( sum( (P_line * ea) .* ea ) );
if normea == 0
    normea = delta;
end

% backstress
a = sum(alphai,2);
norma = sqrt( sum( (P_line * a) .* a ) );
if norma == 0
    norma = delta;
end
% aktueller pseudo Spannungsdeviator
es = P * esig;
% pseudo effektive spannung
ebeta = es  - ea;

% normale
en = P_hat * (ebeta./(sqrt(2)*ek));
transen = en';

%--------------------------------------------------------------------------
%               Ableitung der Teilbackstresstensoren
%--------------------------------------------------------------------------

% normen der teilbackstresstensoren
norm_ai = sqrt( sum( ( P_line * alphai ) .* alphai ) );
norm_ai(norm_ai == 0) = delta;

% normen der pseudo teilbackstresstensoren
norm_eai = sqrt( sum( ( P_line * ealphai ) .* ealphai ) );
norm_eai(norm_eai == 0) = delta;

% hilfsgrößen
Li = alphai./norm_ai;
eLi = ealphai./norm_eai;

% Materialparameter 
c_i = c_i0 .* (ones(1,M) + a_i1.*exp(-b_i1.*p) + a_i2.*exp(-b_i2.*p));
ec_i = ec_i0 .* (ones(1,eM) + ea_i1.*exp(-eb_i1.*p) + ea_i2.*exp(-eb_i2.*p));

chi_i = Q_i .* (2 - transen * P_check * Li) .* (1 + a_chi * exp(b_chi * rm));
echi_i = eQ_i .* (2 - transen * P_check * eLi) .* (1 + ea_chi * exp(eb_chi * erm));

% weitere hilfsgrößen
var1 = c_i .* r_i;
evar1 = ec_i .* er_i;
var2 = A * en;
var3 = norm_ai./r_i;
var3(var3>1) = 1;
evar3 = norm_eai./er_i;
evar3(evar3>1) = 1;

% init ableitungen
dealphai_dp = zeros(ntens,M);
dalphai_dp = zeros(ntens,M);

% Schleife über alle Backstresstensoren
% !!!! BERECHNUNGEN WÜRDE AUCH OHNE SCHLEIFE GEHEN ABER OHNO WANG HAT
% GEZEIGT; DASS ES MIT SCHLEIFE BISSI SCHNELLER IS !!!!
for ii = 1 : M
    
    % backstress
    dalphai_dp(:,ii) = var1(ii) * (var2 - (var3(ii)).^(chi_i(ii)+1).*Li(:,ii));
    
end

for ii = 1 : eM
       
    % pseudo backstress
    dealphai_dp(:,ii) = evar1(ii) * ( var2 - (evar3(ii)).^(echi_i(ii)+1).*eLi(:,ii) );
    
end

da_dp = sum(dalphai_dp,2);
dea_dp = sum(dealphai_dp,2);

%--------------------------------------------------------------------------
%                  Ableitungen der Gedächtniss- und Fließflächen                     
%--------------------------------------------------------------------------

% Hilfsgrößen
L = a./norma;
eL = ea./normea;

dummy = sign(norma - rm);
hg = 0.5 * (dummy + abs(dummy));

dummy = sign(normea - erm);
ehg = 0.5 * (dummy + abs(dummy));

dummy = sum( (P_line * da_dp) .* L);
mla = 0.5 * ( dummy + abs(dummy) );

dummy = sum( (P_line * dea_dp) .* eL);
emla = 0.5 * ( dummy + abs(dummy) );

dummy = 1 - norma/rm;
dummy = 0.5 * ( dummy + abs(dummy) );

edummy = 1 - normea/erm;
edummy = 0.5 * ( edummy + abs(edummy) );

% Ableitung Gedächtnissflächen
drm_dp = hg * mla - cm * dummy;
derm_dp = ehg * emla - ecm * edummy;

% Ableitungen des radius
dek_dp = ek1 * eak * eck * exp(eck*erm) * derm_dp;
% Abfangen von 0*Inf = NaN
dek_dp(isnan(dek_dp)) = 0; 

%--------------------------------------------------------------------------
%                  plastischer tangentenmodul                   
%--------------------------------------------------------------------------

eh = transen * P_check * dea_dp + sqrt(2) * dek_dp;

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
drm = drm_dp .* dp;
dealphai = dealphai_dp .* dp;
derm = derm_dp .* dp;

%--------------------------------------------------------------------------
%                   Zusammenfassen der Inkremente                        
%-------------------------------------------------------------------------- 
  
dX = [deeps; depsp; reshape(dealphai,3*eM,1); dp; derm; reshape(dalphai,3*M,1); drm];

end

