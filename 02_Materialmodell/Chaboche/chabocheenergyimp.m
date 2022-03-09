function [ZVARneu] = chabocheenergyimp(domega,ZVAR,para,dESIG,REFSIG,REFEPS,REFEPSP)
% Chaboche Plastizit‰tsmodell .
% implementiert f¸r die energiegesteuerte integration bei inkrementellen
% Kerbn‰herungen mit euler implizit
%
% INPUT:
%  domega -> Inkrement der energie
%  ZVAR   -> Zustandsvariablen
%  para   -> Parameter des Pseudo Modells und des Materialmodells
%  CEL    -> elastische Steifigkeit
% dESIG   -> pseudo elastisches Inkrement
%  REF... -> Referenzzust‰nde
%
% OUTPUT:
%  ZVARneu -> neuer zustand nach inkrement
%
%__________________________________________________________________________
%
% Zustandsvariablen:
% Zuerst so wie bei dehnungsgesteuerter integration 
% ZVAR = [sig;epsp;alphai;r;p] !! isotrope verfestigung nicht
% ber¸cksichtigt!!
%
% Darstellung von Tensoren
%         sig11              eps11
%  sig =  sig22       eps =  eps22
%         sig12             2eps12
%__________________________________________________________________________          
%
% Parameter:
% zuerst Material- 
% !! isotrope verfestigung nicht ber¸cksichtigt !!!
% ber¸cksichtigt!!
%     para [E,nu,q,gamma,zeta_i,r_i,r0]
%       M = (length(para) - 5)/2
% q = gamma = 0
%
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft                                                   |
% |  Stand: Juli 2020                                                 |
%  ----------------------------------------------------------------------

% -------------------------------------------------------------------------
% !!!!! AKTUELL NUR ESZ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% F¸r weitere Spannungszust‰nde m¸ssten "nur" die Modellgleichung
% implementiert werden
% -------------------------------------------------------------------------
ntens = 3;
ndi = 2;

%--------------------------------------------------------------------------
%                  Materialparamter ermitteln                                  
%--------------------------------------------------------------------------

% Anzahl Backstress
M = (length(para)-5)/2;
% Elastizit‰tskonstanten
E = para(1);
nu = para(2);
r0 = para(end);                                                       % startradius fliessfl‰che
q = 0;
gamma = 0;

% kinematische Verfestigung
c_i = para(5:4+M);
r_i = para(5+M:end-1);
r0 = para(end);                                                       % startradius fliessfl‰che
h_i = c_i.*r_i;

% oft verwendete Konstanden
w3d2 = sqrt(3/2);
w2d3 = sqrt(2/3);

%--------------------------------------------------------------------------
%                        Hilfsmatrizen                                    
%--------------------------------------------------------------------------
% Abbildungen 
[M_, MLINE, MHAT, A, MCHECK] = set_maps(ntens,ndi);                                                                           
DEL = elast_nachgiebigkeit(E,nu,ntens,ndi);                                % Nachgiebigkeit

                          
%--------------------------------------------------------------------------
%                  Zustandsvariablen auslesen                                 
%--------------------------------------------------------------------------

% Spannungen und Dehnungen
sig0 = ZVAR(1:ntens);
epsp0 = ZVAR(ntens+1:2*ntens);
% Backstresstensoren
alpha0 = reshape( ZVAR(2*ntens+1:(M+2)*ntens) , ntens, M);
% plastische bogenl‰nge
p0 = ZVAR(end);

%--------------------------------------------------------------------------
%                  Init Newton iteration                                 
%--------------------------------------------------------------------------
dsig = dESIG;
dp = 0;

%--------------------------------------------------------------------------
%                  kriterien f¸r newton iteration                            
%--------------------------------------------------------------------------
maxiter = 10000;
alpha = 0.0;
tolf = 1e-10;
tole = 1e-16;
iter = 1;



%--------------------------------------------------------------------------
%                  Ausf¸hren erster schritt                              
%--------------------------------------------------------------------------

[s,theta_i,beta,n,depsp,alpha,Fn,R2,normR2] = firstupdate( ...
     domega,M_,sig0, dsig, c_i, dp, alpha0, r0, h_i, MLINE, ...
     MHAT,w3d2,w2d3,A,DEL,REFSIG,REFEPS,epsp0);

%--------------------------------------------------------------------------
%                  newton iter                             
%--------------------------------------------------------------------------
% Iterationsschleife
while Fn > tolf || normR2 > tole
    
    fprintf('F = %.5d R2 = %.5d dp = %.5d dsig = [%.3f %.3f %.3f] dESIG = [%.3f %.3f %.3f] depsp = [%.3d %.3d %.3d]\n',Fn,normR2, dp, dsig, dESIG, depsp);

    % keine endlosschleifen
    iter = iter + 1;
    if iter > maxiter + 1
        msg = ['Keine Konvergenz in Newton/ Verfahren, ',...
               ' Anzahl Iterationen: ', num2str(iter-1),...
               ' Lastschritt: ',num2str(lastschritt),...
               ' Fehlerquadrate:', num2str(norm2)];
		warning(msg)
        break
    end
    
    % berechne Ableitung
    dR1ddp = w3d2 * n' * MCHECK * sum(theta_i.^2 .* c_i.* alpha0,2) ...
            - sum(theta_i.*h_i) ...
            + dp * sum((theta_i.*c_i).^2 .*r_i);
    
    dR1ddsig = w3d2 * n;
    
    dR2ddp = -w3d2 * (sig0 + dsig - REFSIG) .* n - w3d2 * dsig .* n;
    
    dR2ddsig = -diag(sig0 - REFSIG)*DEL - 3/(2*r0) * dp * diag(sig0-REFSIG)*M_ ...
               - 2 * ( diag(DEL*dsig) + diag(dsig)*DEL ) ...
               - 2 * (diag(w3d2*dp*n) + 3/(2*r0) * dp * diag(dsig)*M_) ...
               - diag( DEL*sig0 + epsp0 - REFEPS);
    
    % Berechne neues Inkrement
    dR = [dR1ddp, dR1ddsig'; ...
          dR2ddp, dR2ddsig];
      
    dx = - dR\[Fn;R2];
    
    % update dsig und dp
    dp = dp + dx(1);
    dsig = dsig + dx(2:4);
    
    % update anderes zeug
    [s,theta_i,beta,n,depsp,alpha,Fn,R2,normR2] = firstupdate( ...
     domega,M_,sig0, dsig, c_i, dp, alpha0, r0, h_i, MLINE, ...
     MHAT,w3d2,w2d3,A,DEL,REFSIG,REFEPS,epsp0);
    
end
fprintf('F = %.5d R2 = %.5d dp = %.5d dsig = [%.3f %.3f %.3f] dESIG = [%.3f %.3f %.3f] depsp = [%.3d %.3d %.3d]\n',Fn,normR2, dp, dsig, dESIG, depsp)
fprintf('\n\n\n');
%--------------------------------------------------------------------------
%                  Update Zustandsvariablen                            
%--------------------------------------------------------------------------
% ZVAR = [sig;epsp;alphai;r;p]
% neue variablen
sig = sig0 + dsig;
epsp = epsp0 + depsp;
p = p0 + dp;

% zusammenfassen
ZVARneu = [sig;epsp;alpha(:);r0;p];
end % Ende Funktion 


function [s,theta_i,beta,n,depsp,alpha,Fn,R2,normR2] = firstupdate(domega,M_,sig0, dsig, c_i, dp, alpha0, r0, h_i, MLINE, MHAT,w3d2,w2d3,A,DEL,REFSIG,REFEPS,epsp0)
% Hilfsfunktion um updates zu berechnen
% Spandev
s = M_ * (sig0 + dsig);
% normale
theta_i = 1./(1+c_i.*dp);
beta = r0 * ( s - sum(theta_i .* alpha0,2)) /(r0 + dp * sum(theta_i .* h_i));
nbeta = sqrt( sum( (MLINE*beta) .* beta) );
n = MHAT * beta./nbeta;
% plast dehn inc
depsp = w3d2 * dp * n;
% backstress

alpha = theta_i .* (alpha0 + w2d3 * dp * h_i .* (A * n));
% Flieﬂfl‰che
beta = s - sum(alpha,2);
Fn = w3d2 * sqrt( beta' * MLINE * beta) - r0; 
% Fehler in Neuber
deps = DEL * dsig + depsp;
eps = DEL * ( sig0 + dsig) + epsp0 + depsp;
R2 = domega - (sig0 + dsig - REFSIG) .* deps ...
            - (eps - REFEPS) .* dsig;
normR2 = sum(R2.*R2);
end
