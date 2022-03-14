function [tout,yout] = rk87(odefun, tspan, y0, options, varargin)
% Funktion zur Integration eines Materialmodells mit kinemaitscher und
% isotroper Verfestigung mit explizitem eingebettetem Runkge-Kutta
% -Verfahren achter Ordnung nach:
% P.J. Prince & J.R. Dorman (1981) 
% High order embedded Runge-Kutta formulae. J.Comp. Appl. Math., Vol. 7. p.67-75.
%
% Vorgehen bei Schrittweitensteurung aus Vorlesungsskript:
% Numerischer Methoden der technischen Dynamik, Prof. Schweizer, TU Darmstadt
%
% Weitere Infos zu eingebetten RK-verfahren aus: 
% Numerik gewöhnlicher Differentialgleichungen mit Computer Algebra
% Systemen
%
% Modell muss DGL system der Form y' = f(x,y) sein
% Rückgabewert des Modells (y') muss Spaltenvektor sein
%
%
% INPUT:
%       odefun  -> System gewöhnlicher DGL 1. Ordnung
%       tspan   -> Zeitspanne der Integration ( für ratenunabhängige
%                  Plastizität immer [0,1] )
%       y0      -> Anfangsbedingungen/ momentaner Zustand
%       options -> Optionen (siehe MATLAB Dokumentation)
%                  Verfügbar sind RelTol, AbsTol, MaxStep, InitialStep
%                  wenn keine optionen definiert sind aber zusätzliche
%                  Variablen übergeben werden setze options = []
%       varargin-> Optional , alles hinter options wird DGL übergeben
%                   gedacht für Materialparameter und Lastinkremente
%
% OUTPUT:
%       tout -> Zeitpunkte als Spaltenvektor (nur damit auch allg. DGL
%               integriert werden können, für ratenunabhängige Plastizität
%               nicht nötig)
%       yout -> num. Lösungen der DGl 
%               alle Spalten in einer zeile bilden die Lösung zum
%               Zeitpunkt gespeichert in der Zeile von tout 
%
% ------------------------------------------------------------------------
%
% ------------------------------------------------------------------------
% | Jan Kraft Tu Darmstadt FG Werkstoffmechanik                          |
% | Stand : Juni 2019                                                    |
% ------------------------------------------------------------------------

% festes (fiktives) Intervall t in [0,1] und Initialisieren von t
t0 = tspan(1); % Start des Intervalls
tend = tspan(2); % Ende des Intervalls
t = t0;
hmin = 1e-13; % Festlegen minimaler Schrittweite

% Koef. des Verfahrens
c_i=  [ 1/18, 1/12, 1/8, 5/16, 3/8, 59/400, 93/200, 5490023248/9719169821, 13/20, 1201146811/1299019798, 1, 1]';
a_ij = [  1/18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
          1/48, 1/16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
          1/32, 0, 3/32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
          5/16, 0, -75/64, 75/64, 0, 0, 0, 0, 0, 0, 0, 0, 0; 
          3/80, 0, 0, 3/16, 3/20, 0, 0, 0, 0, 0, 0, 0, 0; 
          29443841/614563906, 0, 0, 77736538/692538347, -28693883/1125000000, 23124283/1800000000, 0, 0, 0, 0, 0, 0, 0;
          16016141/946692911, 0, 0, 61564180/158732637, 22789713/633445777, 545815736/2771057229, -180193667/1043307555, 0, 0, 0, 0, 0, 0;
          39632708/573591083, 0, 0, -433636366/683701615, -421739975/2616292301, 100302831/723423059, 790204164/839813087, 800635310/3783071287, 0, 0, 0, 0, 0;
          246121993/1340847787, 0, 0, -37695042795/15268766246, -309121744/1061227803, -12992083/490766935, 6005943493/2108947869, 393006217/1396673457, 123872331/1001029789, 0, 0, 0, 0;
         -1028468189/846180014, 0, 0, 8478235783/508512852, 1311729495/1432422823, -10304129995/1701304382, -48777925059/3047939560, 15336726248/1032824649, -45442868181/3398467696, 3065993473/597172653, 0, 0, 0;
          185892177/718116043, 0, 0, -3185094517/667107341, -477755414/1098053517, -703635378/230739211, 5731566787/1027545527, 5232866602/850066563, -4093664535/808688257, 3962137247/1805957418, 65686358/487910083, 0, 0;
          403863854/491063109, 0, 0, -5068492393/434740067, -411421997/543043805, 652783627/914296604, 11173962825/925320556, -13158990841/6184727034, 3936647629/1978049680, -160528059/685178525, 248638103/1413531060, 0, 0]';
 b_8 = [ 14005451/335480064, 0, 0, 0, 0, -59238493/1068277825, 181606767/758867731,   561292985/797845732,   -1041891430/1371343529,  760417239/1151165299, 118820643/751138087, -528747749/2220607170,  1/4]';
 b_7 = [ 13451932/455176623, 0, 0, 0, 0, -808719846/976000145, 1757004468/5645159321, 656045339/265891186,   -3867574721/1518517206,   465885868/322736535,  53011238/667516719,                  2/45,    0]';

 % Überprüfe Eingangsdaten
 if nargin < 4
     options = [];
     varargin = {}; % Für Matmodell überflüssig nur für TestDGL
     if nargin < 3
        msg = 'Zu wenig Übergabeparameter in Runge-Kutta';
        error(msg)
     end
 elseif nargin < 5
     varargin = {}; % Für Matmodell überflüssig nur für TestDGL
 end
     
 % Maximale Schrittweite
 hmax = odeget(options, 'MaxStep');
 if isempty(hmax)
     hmax = tend - t0;
 elseif hmax > tend - t0
    hmax = tend - t0;
    msg = ['Maximal zulässige Schrittweite ist' ,...
           ' hmax = 1. Maximaleschrittweite wurde zurückgesetzt!'];
    warning(msg)
 end
 
 % Initiale Schrittweite
 h = odeget(options, 'InitialStep');
 if isempty(h)
%      h = (tend - t0);
     h = (tend - t0)/2; % Versuch mit halber maximalen Schrittweite 
 elseif h > hmax
     h = (tend - t0)/2;
     msg = ['Maximal zulässige Schrittweite ist' ,...
           ' hmax = 1. Schrittweite wurde auf (tend - t0)/2 ',...
           'zurückgesetzt!'];
    warning(msg)
 end
 
 % Realtive Fehlertoleranz (für alle Komponenten des Lösungsvektors)
 rtol = odeget(options, 'RelTol'); 
 if isempty(rtol)
     rtol = 1e-8; % festlegen relativer fehlertoleranz   
 end
 
 % Absolute Fehlertoleranz (für alle Komponenten des Lösungsvektors)
 atol = odeget(options, 'AbsTol'); 
 if isempty(atol)
     atol = 0;%rtol; % festlegen relativer fehlertoleranz                
 end
 
 % Initialisierung
 n_steps = 0; % Schrittzähler
 n_reject = 0; % zahler abgelehnte Schritte
 y = y0(:);    % Startpunkt, abfangen flasches ausgabeformat
 f = y * zeros(1,13); % Vorbelegter Speicher für Funktionsauswertungen
 tout = t;
 yout = y';
 
 % Konstanten für Ablehnung eines Integrationsschrittes
 delta = 0.9; % gewählter Sicherheitsfaktor
 pp1 = 8;     % ordnung des genaueren Verfahrens
 htol = hmin + 1e-13; % Schreittweitentoleranz (nur kontrolle)
  
 % Hauptschleife
 while t < tend && h >= hmin    % schleife wird solange durchlaufen bis Ende
                                % des Integrationsschittes erreicht ist
                                % oder Tolerenazen nicht mehr mit
                                % vernünftiger Schrittweite eingehalten
                                % werden können
      
      % Begrenzen der Schrittweite h wenn Ende des Intervalls erreicht ist                          
      if t + h > tend
          h = tend - t;
          if h < hmin
              disp('kein Fehler letzter schritt nur zu nah an tend, diesen fall noch abfangen')
          end
      end
      % Funktionsauswertungen ( berechne ki des Runge-Kutta-Verfahrens)
      f(:,1) = feval(odefun, t, y, varargin{:});
      for j = 1: 12
          f(:,j+1) = feval(odefun, t+c_i(j)*h, y+h*f*a_ij(:,j), varargin{:});
      end
      % berechnen der zwei Lösungen
      sol_8=y+h*f*b_8;
      sol_7=y+h*f*b_7; 
      % Lokaler Fehlerschätzer
      lok_error = sol_7 - sol_8;
      % Fehler Norm
      err = norm(lok_error,'inf');
      % Toleranz
      tol = atol + rtol * max(norm(sol_7,'inf'),1);
      % neuer Vorschlag für schrittweite
      if err == 0 % teilen durch null abfangen
          err = eps * 100;
      end
      hopt = h * (tol/err)^(1/pp1);
      % Schritt wird akzeptiert wenn err < tol
      % Schrittweitensteuerung
      if err <= tol   
          t = t + h;
          y = sol_8;
          tout = [tout;t];
          yout = [yout;y']; % speichern aller ergebnisse
          %tout = t;
          %yout = y';      % Nur ausgeben des letzten Ergebnisses 
                            % (schneller)
          h = max( min(delta * hopt, hmax), hmin);
      elseif err > tol && h > hmin % schritt wird Verworfen
          n_reject = n_reject + 1;
          h = max(delta * hopt, hmin);
      else % Schritt wurde mit h_min gerechnet und Fehlertoleranz konnte
           % nicht eingehalten werden
          msg = ['angegebene Fehlertoleranzen konnten mit minimaler ' ,...
                 'Schritweite hmin = ', num2str(hmin), ...
                 ' nicht eingehalten werden.'];
          warning(msg);
          return;
      end
      % Kontrolle der Schritweite
      if h - hmin < htol
          fprintf('h = %.3d\n',h)
          fprintf('err = %.3d\n',err)
          msg = ['Schrittweite bei t = ' , num2str(t), ...
                 ' sehr nahe an minimaler Schrittweite'];
          warning(msg)
      end
      n_steps = n_steps + 1;
  end
%  fprintf('Schleifendurchläufe : %i\n',n_steps)
%  fprintf('Abgelehnte Schritte : %i\n',n_reject)
end