function [ZVOUT, dOOUT] = bisektion(ZVIN, dOIN,...                    
                                           REFSIG, REFEPS,...                
                                           DEL, para, matfun,...             
                                           verfahren,...                     
                                           maxiter,tol,alpha,...             
                                           dEEPS,...
                                           lastschritt,...
                                           varargin)
% BisektionvVerfahren für inkrementelle Kerbnäherung. Zuerst Newton Raphson
% Verfahren bis überschießen in einer komponente festegestellt wird. Dann
% Bisektionsverfahren, ohne dass ein weitere Newton Schritt ausgeführt
% wird.
%
% INPUT:
% ZVIN        -> aktuelle Zustandsvariablen e R^(xx,1)
%                ersten beiden Tensoren müssen Spannungen und plastische
%                Dehungen sein
% dOIN        -> Inkrement der Energie e R^(ntens,1)
% REF...      -> Referenzzustände e R^(ntens,1)
% DEL         -> elastische Nachgiebigkeit e R^(ntens x ntens)
% para        -> Modell parameter
% matfun      -> Function handle für Materialfunktion
% verfahren   -> Name des Verfahren, steuert energie- und
%                ableitungsrechnung (str)
% maxiter     -> maximal erlaubte anzahl an iterationen
% tol         -> toleranz für abbruchbedingung
% alpha       -> relaxationskoef.
% dEEPS       -> erstes Inkrement der Dehnungen e R^(ntens,1)
% lastschritt -> aktueller Lastschritt (für fehlerausgabe) e int
% varargin    -> variabler Input für Energie und Ableitungsrechnung
%
%
% OUTPUT:
% ZVOUT       -> neue Zustandsvariablen
% dOOUT       -> Energieinkrement mit neuen Zustandsvariablen
% _________________________________________________________________________


% Spannungszustand (aktuell nur für ESZ)
ntens = 3;
ndi = 2;

% Funktionen je nach verfahren
switch verfahren
    
    case 'Neuber'
        energyfun = @neuberenergy;
        ablfun = @neuberableitung;
    case 'ModNeuber'
        energyfun = @modneuberenergy;
        ablfun = @modneuberableitung;
    case 'ESED'
        energyfun = @esedenergy;
        ablfun = @esedableitung;
    case 'UniExp'
        energyfun = @uniexpenergy;
        ablfun = @uniexpableitung;
end


% Init iteration mit elastischem Dehnungsinkrement
dEPS = dEEPS;
[ZVOUT,DEP,CEP] = matfun(ntens,ndi,dEPS,ZVIN,1,para);

% Rauslesen der Zustände
SIG = ZVOUT(1:ntens);
EPS = DEL * SIG + ZVOUT(1+ntens:2*ntens);
dEPS = EPS - DEL * ZVIN(1:ntens) - ZVIN(ntens+1:2*ntens);
dSIG = SIG - ZVIN(1:ntens);

% berechne Fehler
dOOUT = energyfun(SIG,dSIG,EPS,dEPS,REFSIG,REFEPS,varargin{:});
err = dOIN - dOOUT;
norm2 = sum(err.*err);

% Iterationsschleife
fiter = zeros(3,maxiter+1);           % speicher zum debuggen
depsiter = zeros(3,maxiter+1);        % Speicher zum debuggen
dsigiter = zeros(3,maxiter+1);        % Speicher zum debuggen
FP = zeros(3,1);                      % Zielfunktion an Plus Grenze
FM = zeros(3,1);                      % Zielfunktion an Minus Grenze
DEP = zeros(3,1);                     % Inkrement an Plus Grenze
DEM = zeros(3,1);                     % Inkrement an Minus Grenze
tol_overshot = 1e-1;                   % Toleranz grenze fuer ueberschießen
verfahren = [ 0, 0, 0 ];              % verwendetes iterationsverfahren 
                                      % 0 = Newton
                                      % 1 = Bisektion
iter = 1;

while norm2 > tol
        
    % keine endlosschleifen
    iter = iter + 1;
    if iter > maxiter + 1
        msg = ['Keine Konvergenz in Bisektions Verfahren, ',...
               ' Anzahl Iterationen: ', num2str(iter-1),...
               ' Lastschritt: ',num2str(lastschritt),...
               ' Fehlerquadrate:', num2str(norm2)];
		warning(msg)
        break
    end
    
    % Speicher zum Debuggen
    depsiter(:,iter) = dEPS;
    dsigiter(:,iter) = dSIG;
    fiter(:,iter) = err;
    
    if sum(verfahren) < 2 % In Mindestens einem Verfahren wird noch Newton 
                          % Verfahren angewandt.
        % Ableitung des Näherungsvrefahrens
        [G,P] = ablfun(SIG,EPS,REFSIG,REFEPS,varargin{:});
        derr = G * CEP + P;

        % Änderung des Dehnungsinkrements
        ddEPS = derr\err;

        % dämpfen der Schrittweite
        s = 1 - alpha;         % konstante Dämpfung

        % neues (relaxiertes) Dehnungsinkrement nach Newton
        dEPSnewton = dEPS + s .* ddEPS;
    end
    
    % Neues Dehnungsinkrement nach newton oder bisektion
    for ee = 1 : ntens
        if verfahren(ee) % Bisektion
            
            % Neues Dehungsinkrement
            [dEPS(ee),FP(ee),FM(ee),DEP(ee),DEM(ee)] = bisek(...
                             FP(ee),FM(ee),...
                             DEP(ee),DEM(ee),...
                             fiter(ee,iter),dEPS(ee)); 
        else % Newton
            
            % teste auf ueberschiesen
            [dEPS(ee),verfahren(ee),FP(ee),FM(ee),DEP(ee),DEM(ee)] = ...
                overshoot(fiter(ee,iter-1),fiter(ee,iter),...
                          depsiter(ee,iter-1),depsiter(ee,iter),...
                          dEPSnewton(ee),...
                          tol_overshot);
        end
    end
     
    % Integration mit neuem Dehnungsinkrement
    [ZVOUT,DEP,CEP] = matfun(ntens,ndi,dEPS,ZVIN,1,para);
    
    % Rauslesen der Zustände
    SIG = ZVOUT(1:ntens);
    EPS = DEL * SIG + ZVOUT(1+ntens:2*ntens);
    dSIG = SIG - ZVIN(1:ntens);

    % berechne Fehler
    dOOUT = energyfun(SIG,dSIG,EPS,dEPS,REFSIG,REFEPS,varargin{:});
    err = dOIN - dOOUT;
    norm2 = sum(err.*err);

end

% Speicher zum Debuggen
depsiter(:,iter+1) = dEPS;
dsigiter(:,iter+1) = dSIG;
fiter(:,iter+1) = err;


end % Ende Hauptfunktion

% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %
% ::           HILFSFUNKTION Bisektion                                 :: %
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %

function [deps,fp,fm,dep,dem] = bisek(fp,fm,dep,dem,f,deps) 
    % Prüft ob Bisektion ausgeführt werden muss und führt sie
    % gegebenenfalls aus
    % INPUT:
    % fp       -> zielfun plus grenze
    % fm       -> zielfun minus grenz
    % dep      -> Inkrement plus grenze
    % dem      -> inkrement minus grenze
    % f        -> zielfun aktuell
    % deps     -> Inkrement aktuell
    
    % OUTPUT:
    % deps -> neues/oder alters Dehnungsinkrement
    % .... -> neue grenzen
    
    % prüfe ob f zielfunktion verbessert
%     if abs(f) > min(abs(fp),abs(fm))
%         msg = 'Bisektion gescheitert';
%         error(msg)
%     end
    
    % f > 0
    if f > 0 
        
        % neue obere grenzen
        fp = f;
        dep = deps;
        
    else
        
        % neue untere grenzen
        fm = f;
        dem = deps;
        
    end
    
    % neues Dehnungsinkrement
    dist = -fp/(fm-fp);
    deps = dist * dem + (1-dist) * dep;
        


end % Ende Bisektion


% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %
% ::           HILFSFUNKTION überschießen                              :: %
% ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: %
function [deps,verfahren,fp,fm,dep,dem] = overshoot( ...
                      f1,f2,deps1,deps2,depsnewton,tol ) 
    % Prüft ob Bisektion ausgeführt werden muss und führt sie
    % gegebenenfalls aus
    % INPUT:
    % f1         -> alte zielfunktion in iterationschritt: iter - 1
    % f2         -> neue zielfunktion in iterationschritt: iter
    % deps1      -> altes Inkrement in iterationschritt  : iter - 1
    % deps2      -> neues Inkrement in iterationschritt  : iter
    % depsnewton -> dehnungsinkrement aus newton schritt
    % tol        -> toleranz welche Potenzunterschiede als ueberschießen
    %               bewertet werden
    % OUTPUT:
    % deps       -> neues/oder alters Dehnungsinkrement
    % verfahren  -> 0 newton, 1 bisektion
    % fp,fm      -> zielfunktion plus und minus
    % dep,dem    -> Inkremente plus und minus
    
    % prüfe vorzeichenwechsel in zielfunktion
    if f1 * f2 < 0
        
        % Prüfe ob vorzeichenwechsel außerhalb der Tolranz ist (VZ wechsel
        % aber trozdem Konvergenz)
        if abs(f2/f1) > tol
            
           % Anstatt Bisektion, Interpolation linear fuer bessere
           % Konvergenz
           dis = -f1 / (f2-f1);
           
           % neues Dehnungsinkrement
           deps = dis * deps2 + (1-dis) * deps1;
           
           % Verfahren umstellen auf bisektion
           verfahren = 1;
           
           % feststellen der Plus und Minus Grenzen
           if f1 > 0 % f1 ist plus grenze
               
               fp = f1;
               fm = f2;
               dep = deps1;
               dem = deps2;
               
           else     % f2 ist plus grenze
               
               fm = f1;
               fp = f2;
               dem = deps1;
               dep = deps2;
               
           end
        else % VZ wechsel innerhalb toleranz -> newton schritt wird angenommen
            
            deps = depsnewton;
            verfahren = 0;
            fp = 0;
            fm = 0;
            dep = 0;
            dem = 0;
            
        end
        
    else % kein vz wechsel -> newton schritt wird angenommen
        
        deps = depsnewton;
        verfahren = 0;
        fp = 0;
        fm = 0;
        dep = 0;
        dem = 0;
        
    end

end % Ende Overschoot