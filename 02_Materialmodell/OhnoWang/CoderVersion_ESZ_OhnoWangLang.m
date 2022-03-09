%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Matlab Coder Version von:                                            %
%    Modellgleichung ESZ Ohno Wang Modell für Integration nach Lang       %
%                                                                         %
%    Aufgerufen in:                                                       %
%    ohnowanglang.m                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dX = CoderVersion_ESZ_OhnoWangLang(~, X, M, eM, desig, para, epara, CEL, DEL, ...
                        P, P_line, P_hat, A, P_check)
% Inkrements des System von DGL fürs Ohno Wang Modell mit den Ansatz von Lang
%
% INPUT:
%      ~       -> Zeit wird hier nicht gebraucht, aber wird von ode solver
%                 benötigt
%     X      -> Speicher für aktuellen Zustand, begrenzt auf 20
%                 Backstresstensoren 
%  para       -> Speicher für Material- und Struckturparameter
%                 begrenzt auf 20 Backstresstensoren
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
Li = zeros(3,M);
eLi = zeros(3,eM);
var3 = zeros(1,M);
evar3 = zeros(1,eM);
for kk = 1 : M
    Li(:,kk) = alphai(:,kk)/norm_ai(kk);
    var3(kk) = norm_ai(kk)/r_i(kk);
end

for kk = 1 : eM
    eLi(:,kk) = ealphai(:,kk)/norm_eai(kk);
    evar3(kk) = norm_eai(kk)/er_i(kk);
end

var1 = c_i .* r_i;
evar1 = ec_i .* er_i;
var2 = A * en;
var3(var3>1) = 1;
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
