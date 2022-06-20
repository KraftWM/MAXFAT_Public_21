function [SIG, EPS, EPSP, ALPHA, R, P] = esedimp( ...
                                   ZVAR0, para, ESIG, EEPS, ...
                                   ntens, ndi, material, numink)
% Implentation der ESED Methode f�r mehrachsige Kerbn�herung
%
%  !!!!!    Implementierung f�r implizite Integration !!!!!!!
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
        matfun = @chaboche;
    case 'OhnoWang'
        matfun = @ohnowang;
    case 'Jiang'
        matfun = @jiang;
    otherwise
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)
end 

% -------------------------------------------------------------------------
% Elastizit�tskonstanten
E = para(1);                                                               % E-Modul
nu = para(2);                                                              % Querdehnzahl
% CEL = elast_steifigkeit(E,nu,ntens,ndi);                                 % Steifigkeit
DEL = elast_nachgiebigkeit(E,nu,ntens,ndi);                                % Nachgiebigkeit

% -------------------------------------------------------------------------
% Init Referenzzust�nde
% Annhame Lastfolge startet im Nullzustand -> Alle Refs mit 0
% inititalisieren
% REF(3,3) = [REF_ESIG , REF_EEPS, REF_SIG] -> Referenzpunkte in plast.
% Dehnungen werden net gebraucht
REF = zeros(3,3);

% -------------------------------------------------------------------------
% Speicher f�r Zustandsvariablen
ZVAR = zeros(size(ZVAR0,1),numink);
ZVAR(:,1) = ZVAR0;
% ZVAR(1:ntens,1:jj) = ESIG(:,1:jj);

% -------------------------------------------------------------------------
% Init altes Inkrement der pseudo Gr��en
dESIG0 = zeros(3,1);

% -------------------------------------------------------------------------
% konstanten f�r iterationschleife
TOL = 1e-16;        % Toleranz
maxiter = 10;      % Maximale iterationen (keine endlosschleife)


% -------------------------------------------------------------------------
% Hauptschleife �ber alle Inkremente der Lastfolge
for ii = 2 : numink
   
    % =====================================================================
    % Inkremente in Pseudo elastischen Gr��en
    dESIG = ESIG(:,ii) - ESIG(:,ii-1);
    dEEPS = EEPS(:,ii) - EEPS(:,ii-1);
    
    % =====================================================================
    % Pr�fe Umkehrpunkt und updatete ggf. die Referenzzust�nde
    idx = dESIG .* dESIG0 < 0;
    REF(idx,:) = [ESIG(idx,ii-1), EEPS(idx,ii-1), ZVAR(idx,ii-1)];
    
    % =====================================================================
    % Speichern pseudo Spannungsinkrement falls nicht gleich 0
    idx = dESIG ~= 0;
    dESIG0(idx) = dESIG(idx);
    
    % =====================================================================
    % Inkrement der pseudo Verzerrungsenergiedichte an der stelle ii
    % dEESED(i) = ( ESIG(i)+ dESIG - REF_ESIG(i) ) * dEEPS(i)
    dEESED = (ESIG(:,ii-1) - REF(:,1)).* dEEPS + dEEPS .* dESIG;
    
    % =====================================================================
    % erster elastischer Versuch
    dEPS = dEEPS;
    % Integration MatMod mit elastischem Inkrement (stellt akt. Zustand
    % her)
    [ZVAR(:,ii)] = matfun(ntens,ndi,dEPS,ZVAR(:,ii-1),1,para);
    % Inkrement der Verzerrungsenergiedichte an stelle ii
    dSIG = ZVAR(1:ntens,ii) - ZVAR(1:ntens,ii-1);
    dESED = (ZVAR(1:ntens,ii-1) - REF(:,3)) .* dEPS + dEPS .* dSIG;
    % Zielfunktion
    f = dEESED - dESED;
    % Norm
    fnorm2 = sum(f.*f);
    
    % =====================================================================
    % Iterationsschleife
    iter = 0;                 % zahler iterationsschleifen
    store_fnorm = zeros(1,101);  % speicher f�r Norm zum Debuggen
    store_f = zeros(3,101);  % speicher Zielfun zum Debuggen
    store_deps = zeros(3,101);   % speicher f�r deps zum debuggen
    store_ddeps = zeros(3,100);  % speicher f�r �nderung deps
    
    while fnorm2 > TOL
        
        % keine Endlosschleifen
        if iter == maxiter
            msg = ['keine Konvergenz in impliziter ESED Methode',...
                   ' in Lastschritt ', num2str(ii)];
            error(msg)
        end
        
        % Matrix der Ableitungen an der Stelle ii
        Pii = diag(ZVAR(1:ntens,ii) - REF(:,3));
        dfdde = -Pii;
        % Abfangen singul�re Matrix
        for i = 1:ntens
            if dfdde(i,i) == 0
                dfdde(i,i) = 1;
            end
        end
        
        % �nderung in Dehnungsinkrement
        ddeps = -dfdde\f;
        
        % speichern zum debuggen
        store_fnorm(iter+1) = fnorm2;
        store_f(:,iter+1) = f;
        store_deps(:,iter+1) = dEPS;
        store_ddeps(:,iter+1) = ddeps;
        
        % neues Dehnungsinkrement
        dEPS = dEPS + 0.8* ddeps;
        
        % Integration MatMod
        [ZVAR(:,ii)] = matfun(ntens,ndi,dEPS,ZVAR(:,ii-1),1,para);
        
        % Inkrement der Verzerrungsenergiedichte an stelle ii
        dSIG = ZVAR(1:ntens,ii) - ZVAR(1:ntens,ii-1);
        dESED = (ZVAR(1:ntens,ii-1) - REF(:,3)) .* dEPS + dEPS .* dSIG;
        
        % Zielfunktion
        f = dEESED - dESED;
        
        % Norm
        fnorm2 = sum(f.*f);
        
        % Displayausgabe
%         fprintf('iter: %.2d   fnorm: %.2d\n',iter,fnorm);
        
        % Inkrement iterationsz�hler
        iter = iter + 1;

    end % Ende Itertionsschleife
    
end % Ende Schleife �ber alle Inkremente



% -------------------------------------------------------------------------
% Herrauslesen der lokalen Gr��en je nach Material
switch material
    case 'Chaboche'
        
        M = (length(para)-5)/2;  
        SIG = ZVAR(1:ntens,:);
        EPSP = ZVAR(ntens+1:2*ntens,:);
        EPS = EPSP + DEL * SIG;
        ALPHA = zeros(ntens,M,numink);
        for i = 1:M
            ALPHA(:,i,:) = ZVAR((i+1)*ntens+1:(i+2)*ntens,:);
        end
        R = ZVAR(end-1,:);
        P = ZVAR(end,:);
        
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
        
        M = (length(para)-8)/7;  
        SIG = ZVAR(1:ntens,:);
        EPSP = ZVAR(ntens+1:2*ntens,:);
        EPS = EPSP + DEL * SIG;
        ALPHA = zeros(ntens,M,numink);
        for i = 1:M
            ALPHA(:,i,:) = ZVAR((i+1)*ntens+1:(i+2)*ntens,:);
        end
        P = ZVAR(end-1,:);
        R = para(end) * ones(1,numink);
        
    otherwise
        msg = 'Angegebenes Materialmodell nicht implementiert';
        error(msg)
end



end % Ende Funktion
