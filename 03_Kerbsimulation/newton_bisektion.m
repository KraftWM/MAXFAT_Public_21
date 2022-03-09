function [ZVOUT, dOOUT] = newton_bisektion(ZVIN, dOIN,...                    
                                           REFSIG, REFEPS,...                
                                           DEL, para, matfun,...             
                                           verfahren,...                     
                                           maxiter,tol,alpha,...             
                                           dEEPS,...
                                           lastschritt,...
                                           varargin)
% Newton-Raphson Verfahren f�r inkrementelle Kerbn�herung. wenn
% �berschie�en in einer komponente festegestellt wird, wird nur das neue
% Inkrement mit dem Bisektionsverfahren berechnet
%
% INPUT:
% ZVIN        -> aktuelle Zustandsvariablen e R^(xx,1)
%                ersten beiden Tensoren m�ssen Spannungen und plastische
%                Dehungen sein
% dOIN        -> Inkrement der Energie e R^(ntens,1)
% REF...      -> Referenzzust�nde e R^(ntens,1)
% DEL         -> elastische Nachgiebigkeit e R^(ntens x ntens)
% para        -> Modell parameter
% matfun      -> Function handle f�r Materialfunktion
% verfahren   -> Name des Verfahren, steuert energie- und
%                ableitungsrechnung (str)
% maxiter     -> maximal erlaubte anzahl an iterationen
% tol         -> toleranz f�r abbruchbedingung
% alpha       -> relaxationskoef.
% dEEPS       -> erstes Inkrement der Dehnungen e R^(ntens,1)
% lastschritt -> aktueller Lastschritt (f�r fehlerausgabe) e int
% varargin    -> variabler Input f�r Energie und Ableitungsrechnung
%
%
% OUTPUT:
% ZVOUT       -> neue Zustandsvariablen
% dOOUT       -> Energieinkrement mit neuen Zustandsvariablen
% _________________________________________________________________________


% Spannungszustand (aktuell nur f�r ESZ)
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
[ZVOUT,~,CEP] = matfun(ntens,ndi,dEPS,ZVIN,1,para);

% Rauslesen der Zust�nde
SIG = ZVOUT(1:ntens);
EPS = DEL * SIG + ZVOUT(1+ntens:2*ntens);
dEPS = EPS - DEL * ZVIN(1:ntens) - ZVIN(ntens+1:2*ntens);
dSIG = SIG - ZVIN(1:ntens);

% berechne Fehler
dOOUT = energyfun(SIG,dSIG,EPS,dEPS,REFSIG,REFEPS,varargin{:});
err = dOIN - dOOUT;
norm2 = sum(err.*err);

% Speicher zum Debuggen
fiter = zeros(3,maxiter+1);           % speicher zum debuggen
depsiter = zeros(3,maxiter+1);        % Speicher zum debuggen
dsigiter = zeros(3,maxiter+1);        % Speicher zum debuggen
tol_overshot = 1e-1;                   % Toleranz grenze fuer ueberschie�en
iter = 1;

% Iterationsschleife
while norm2 > tol
        
    % keine endlosschleifen
    iter = iter + 1;
    if iter > maxiter + 1
        msg = ['Keine Konvergenz in Newton/Bisektions Verfahren, ',...
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
    
    % Ableitung des N�herungsvrefahrens
    [G,P] = ablfun(SIG,EPS,REFSIG,REFEPS,varargin{:});
    derr = G * CEP + P;
    
    % �nderung des Dehnungsinkrements
	ddEPS = derr\err;
    
    % d�mpfen der Schrittweite
    s = 1 - alpha;         % konstante D�mpfung
    
    % neues (relaxiertes) Dehnungsinkrement nach Newton
	dEPS = dEPS + s .* ddEPS;
    
    % Bisektion bei �berschie�en
    for ee = 1 : ntens
        dEPS(ee) = bisek(fiter(ee,iter-1),fiter(ee,iter),...
                         depsiter(ee,iter-1),depsiter(ee,iter),...
                         dEPS(ee),...
                         tol_overshot); 
    end
     
    % Integration mit neuem Dehnungsinkrement
    [ZVOUT,DEP,CEP] = matfun(ntens,ndi,dEPS,ZVIN,1,para);
    
    % Rauslesen der Zust�nde
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

function deps = bisek(f1,f2,deps1,deps2,depsnewton,tol) 
    % Pr�ft ob Bisektion ausgef�hrt werden muss und f�hrt sie
    % gegebenenfalls aus
    % INPUT:
    % f1         -> alte zielfunktion in iterationschritt: iter - 1
    % f2         -> neue zielfunktion in iterationschritt: iter
    % deps1      -> altes Inkrement in iterationschritt  : iter - 1
    % deps2      -> neues Inkrement in iterationschritt  : iter
    % depsnewton -> dehnungsinkrement aus newton schritt
    % tol        -> toleranz welche Potenzunterschiede als ueberschie�en
    %               bewertet werden
    % OUTPUT:
    % deps -> neues/oder alters Dehnungsinkrement
    
    % pr�fe vorzeichenwechsel in zielfunktion
    if f1 * f2 < 0
        
        % Pr�fe ob vorzeichenwechsel au�erhalb der Tolranz ist (VZ wechsel
        % aber trozdem Konvergenz)
        if abs(f2/f1) > tol
            
           % Anstatt Bisektion, Interpolation linear fuer bessere
           % Konvergenz
           dis = -f1 / (f2-f1);
           
           % neues Dehnungsinkrement
           deps = dis * deps2 + (1-dis) * deps1;
            
        else % VZ wechsel innerhalb toleranz -> newton schritt wird angenommen
            
            deps = depsnewton;
        end
        
    else % kein vz wechsel -> newton schritt wird angenommen
        
        deps = depsnewton;
        
    end


end % Ende Bisektion