function [SIG, EPS, EPSP, ALPHA, R, P] = esedpropeps( ...
                                   ZVAR0, para, ESIG, EEPS, ...
                                   ntens, ndi, material, numink)
% Implentation der ESED Methode f�r mehrachsige Kerbn�herung
%
%  !!!!!    Implementierung sodass dehnungen proportional zu 
%           PseudoDehnungen sind !!!!!!!
%
% INPUT:
% ZVAR0     -> Startwerte der Zustandsvariablen, Zustandsvariblen f�r die
%              Materialmodelle wie bei dehnungsgesteuerter integration
% para      -> Parameter Materialmodel
% ESIG      -> Pseudo elastische Spannungen
% EEPS      -> Pseudo elastische Dehnungen
% ntens     -> Anzahl Tensorkomponenten
% ndi       -> Anzahl Hauptdiagonalelemente
% material  -> Definiert welches Materialmodell verwendet wird
% numink    -> Anzahl inkremente der Lastfolge
%
% OUTPUT:
% SIG    -> elastisch-plastischer Spannungsverlauf
% EPS    -> elastisch-plastischer Dehnungsverlauf
% EPSP   -> plastischer Dehnungsverlauf
% ALPHA  -> Verlauf der Backstresstensoren
% R      -> Radius FF
% P      -> plastische Bogenl�nge
%
% ------------------------------------------------------------------------
% Autor: Jan Kraft                                                        |
% Stand: Januar 2020                                                      |
% ------------------------------------------------------------------------



% -------------------------------------------------------------------------
% Pr�fe ob ebener Spannungszustand

if ntens ~= 3 && ndi ~= 2
    msg = 'ESED akt. nur f�r ESZ gedacht';
    error(msg)
end

% -------------------------------------------------------------------------
% materialfunktion

switch material
    case 'Chaboche'
        
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)
        
    case 'OhnoWang'
        
%         fprintf('noch ohnowang modell implementieren\n');
        matfun = @ohnowang;
        
    case 'Jiang'
        
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)
        
    otherwise
        
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)
end

% -------------------------------------------------------------------------
% Elastizit�tskonstanten
E = para(1);                                                               % E-Modul
nu = para(2);                                                              % Querdehnzahl
% CEL = elast_steifigkeit(E,nu,ntens,ndi);                                   % Steifigkeit
DEL = elast_nachgiebigkeit(E,nu,ntens,ndi);                                % Nachgiebigkeit

% -------------------------------------------------------------------------
% Init Referenzzust�nde
% Annhame Lastfolge startet im Nullzustand -> Alle Refs mit 0
% inititalisieren
% REF(3,3) = [REF_ESIG , REF_EEPS, REF_SIG] -> Referenzpunkte in plast.
% Dehnungen werden net gebraucht
REF = zeros(3,3);

% -------------------------------------------------------------------------
% Finde ersten plastischen Zustand
weiter = 1;
M = [2,-1,0;-1,2,0;0,0,3]./3;
ML = [2,1,0;1,2,0;0,0,2];
jj = 1;
r = para(end);
while weiter
    s = M * ESIG(:,jj);
    F = s'*ML*s - 2/3*r^2;
    if F < 0
        jj = jj + 1;
    else
        jj = jj - 1;
        weiter = 0;
    end
end

% -------------------------------------------------------------------------
% Speicher f�r Zustandsvariablen
ZVAR = zeros(size(ZVAR0,1),numink);
% ZVAR(:,1) = ZVAR0;
ZVAR(1:ntens,1:jj) = ESIG(:,1:jj);

% -------------------------------------------------------------------------
% Init altes Inkrement der pseudo Gr��en
dESIG0 = zeros(3,1);

% -------------------------------------------------------------------------
% Speichern der Energie (zum testen)
% EESED = zeros(3,numink);

% -------------------------------------------------------------------------
% Speicher Differenz aus soll Energieinkrement und ist Inkrement
% f = zeros(3,numink);
% dESED = zeros(3,1);

% -------------------------------------------------------------------------
% speicher f�r proportionalit�tsfaktoren
rho = ones(1,numink);

% -------------------------------------------------------------------------
% Variablen f�r die Interation
maxiter = 100;
tol = 1e-16;

% -------------------------------------------------------------------------
% Hauptschleife �ber alle Inkremente der Lastfolge
for ii = jj + 1 : numink
    
    % =====================================================================
    % Inkremente in Pseudo elastischen Gr��en
    dESIG = ESIG(:,ii) - ESIG(:,ii-1);
    dEEPS = EEPS(:,ii) - EEPS(:,ii-1);
    
    % =====================================================================
    % Pr�fe auf Umkehrpunkt
    ukp = dESIG .* dESIG0 < 0;
    
    % =====================================================================
    % Update die Referenzzust�nde
    REF(ukp,:) = [ESIG(ukp,ii-1), EEPS(ukp,ii-1), ZVAR(ukp,ii-1)];
    
    
    % =====================================================================
    % Pr�fe auf erste Ecke 
    eck = [false,false,false];
    for ee = 1 : 3
        if dESIG(ee) ~= 0 && dESIG0(ee) == 0 && ii > jj + 1
            eck(ee) = true;
        end
    end
    
    % =====================================================================
    % Speichern pseudo Spannungsinkrement falls nicht gleich 0
    idx = dESIG ~= 0;
    dESIG0(idx) = dESIG(idx);
    
    % =====================================================================
    % Inkrement ider pseudo Verzerrungsenergiedichte
    % dEESED(i) = ( ESIG(i) - REF_ESIG(i) ) * dEEPS(i)
    dEESED = (ESIG(:,ii-1) - REF(:,1)).* dEEPS + 0.5 .* dESIG .* dEEPS;
    deesed = sum(dEESED);

     % Speichern der Energie (zu testzwecken)
%     EESED(:,ii) = EESED(:,ii-1) + dEESED;
    
    % =====================================================================
    % Newtoniteration
    rat = rho(ii-1); % strain ratio
    ratiter = zeros(1,maxiter);
    ratiter(1) = rat;
    iter = 0;
    fnorm2 = 1;
    fiter=zeros(1,maxiter);
    verfahren = true; % true = newton, false = Bisektion
    while fnorm2 > tol
        

        % keine endlosschleifen
        iter = iter + 1;
        if iter > maxiter + 1
            msg = ['Keine Konvergenz in Newton Verfahren, Lastschritt ',...
				       num2str(ii), '. Fehlerquadrate:', num2str(fnorm2)];
			warning(msg)
            break
        end
        
        % speichern verh�ltnis
        ratiter(iter) = rat;
        
        % Dehnungsinkrement
        dEPS = rat .* dEEPS;
        
        % integration Dehnungsgesteuert
        [ZVAR(:,ii)] = matfun(3, 2, dEPS, ZVAR(:,ii-1), 1, para);
        
        % Spannungsinkrement
		dSIG = ZVAR(1:ntens, ii) - ZVAR(1:ntens, ii-1);
        
        % Energieinkrement
        dESED = (ZVAR(1:ntens,ii-1) - REF(:,3)).* dEPS  + 0.5 .* dSIG .* dEPS;
        desed = sum(dESED);
        
        % Zielfunktion
        f = deesed - desed;
        fiter(iter) = f;
        
        % Abbruchkriterium
        fnorm2 = f^2;
        
        % Verzweigung Newton/Bisektion
        if fnorm2 > tol && verfahren % newton
            
            % Ableitung Zielfunktion
            dfdrat = sum((ZVAR(1:ntens,ii-1) + dSIG - REF(:,3)).* dEEPS);

            % inkrement verh�ltnis
            drat = f/dfdrat;

    %         % neues Verh�ltnis
    %         rat = rat + drat;
    %         
            % Bisektion falls newton osziliert
            if iter >= 2
                if fiter(iter) * fiter(iter-1) < 0 % Newton Schritt wird abgelehnt
                    dis = -fiter(iter-1) / (fiter(iter)-fiter(iter-1));
                    rat = dis * rat + (1-dis)*ratiter(iter-1);
                else 
                    rat = rat + drat;
                end
            else
                rat = rat + drat;
            end
        
        elseif fnorm2 > tol % Bisektion
            
        end % Ende verzweigung newton/Bisektion

    end % Ende Iterationsschleife
    
    % Speichern verh�ltnis
    rho(ii) = rat;
    
end % Ende Schleife �ber Inkremente 


% -------------------------------------------------------------------------
% Herrauslesen der lokalen Gr��en je nach Material
switch material
    case 'Chaboche'
        
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)        
%         EPS = ZVAR(1:ntens,:);
%         EPSP = ZVAR(ntens+1:2*ntens,:);
%         ALPHA = zeros(ntens,M,numink);
%         for i = 1:M
%             ALPHA(:,i,:) = ZVAR((i+1)*ntens+1:(i+2)*ntens,:);
%         end
%         R = ZVAR(end-1,:);
%         P = ZVAR(end,:);
        
    case 'OhnoWang'
        
        M = (length(para)-3)/3;
        SIG = ZVAR(1:ntens,:);
        EPSP = ZVAR(ntens+1:2*ntens,:);
        EPS = EPSP + DEL * SIG;
        ALPHA = zeros(ntens,M,numink);
        for i = 1:M
            ALPHA(:,i,:) = ZVAR((i+1)*ntens+1:(i+2)*ntens,:);
        end
        P = ZVAR(end,:);
        R = para(end) * ones(1,numink);
        
    case 'Jiang'
        
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)
%         EPS = ZVAR(1:ntens,:);
%         EPSP = ZVAR(ntens+1:2*ntens,:);
%         ALPHA = zeros(ntens,M,numink);
%         for i = 1:M
%             ALPHA(:,i,:) = ZVAR((i+1)*ntens+1:(i+2)*ntens,:);
%         end
%         P = ZVAR(end-1,:);
%         R = para(end) * ones(1,numink);
        
    otherwise
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)
end




end