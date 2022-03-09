function [SIGV , EPSV] = vonmisesvergleich(SIG,EPSE,EPSP,nu,ntens,varargin)
% Funktion berechnet von Mises Vergleichsspannung  
% Spannungszustand, Tensoren in ABAQUS Notation
% INPUT:
%       SIG      -> Spannungen 
%       EPSE     -> elastische Dehnung
%       EPSP     -> plastische Dehnung
%       ntens    -> Anzahl SpanTensoreinträge (Unterscheidung 3D, ESZ, EVZ)
%       varargin -> Variabler input
%
%
% OUTPUT:
%       SIGV,EPSV -> von Mises Vergleichsspannung und -dehnung
%
%
% Notationen:
%  
%           sig_11              epse_11            
%           sig_22              epse_22         
%   sig =   sig_33       epse = epse_33           
%           sig_32              2*epse_32     
%           sig_13              2*epse_13      
%           sig_12              2*epse_12            
%  !!! In Dehnungen wird mit technischen gleitungen gerechnet -> also 
%      doppelten tensorkomponenten !!!
%  ----------------------------------------------------------------------
% |  Autor: Jan Kraft 2019 Tu Darmstadt FG Werkstoffmechanik            |
% |  Stand: Juni 2019                                                   |
%  ----------------------------------------------------------------------



% Fehler abfangen




% Speicher
SIGV = zeros(size(SIG,2),1);
EPSV = SIGV;
% statische Matrizen
P_tilde = diag([1,1,1,0.5,0.5,0.5]);
% ausrechnen Vergleichspannungen
for i = 1 : size(SIG,2)
    
    if ntens == 6 % 3D Spannungszustand
        sig = SIG(:,i);
        epsp = EPSP(:,i);
        epse = EPSE(:,i);
    elseif ntens == 2 % Sigma-Tau
        sig = zeros(6,1);
        epsp = zeros(6,1);
        epse = zeros(6,1);
        sig(1) = SIG(1,i);
        sig(2) = 0;
        sig(3) = 0;
        sig(4) = SIG(2,i);
        sig(5) = 0;
        sig(6) = 0;
        epse(1) = EPSE(1,i);
        epse(2) = -nu*EPSE(1,i);
        epse(3) = -nu*EPSE(1,i);
        epse(4) = EPSE(2,i);
        epse(5) = 0;
        epse(6) = 0;
        epsp(1) = EPSP(1,i);
        epsp(2) = -0.5*EPSP(1,i);
        epsp(3) = -0.5*EPSP(1,i);
        epsp(4) = EPSP(2,i);
        epsp(5) = 0;
        epsp(6) = 0;
    elseif ntens == 3 % ESZ
        sig = zeros(6,1);
        epsp = zeros(6,1);
        epse = zeros(6,1);
        sig(1) = SIG(1,i);
        sig(2) = SIG(2,i);
        sig(3) = 0;
        sig(4) = SIG(3,i);
        sig(5) = 0;
        sig(6) = 0;
        epse(1) = EPSE(1,i);
        epse(2) = EPSE(2,i);
        epse(3) = -nu/(1-nu)*(EPSE(1,i)+EPSE(2,i));
        epse(4) = EPSE(3,i);
        epse(5) = 0;
        epse(6) = 0;
        epsp(1) = EPSP(1,i);
        epsp(2) = EPSP(2,i);
        epsp(3) = -(EPSP(1,i)+EPSP(2,i));
        epsp(4) = EPSP(3,i);
        epsp(5) = 0;
        epsp(6) = 0;
%         msg = 'ESZ noch nicht implementiert';
%         error(msg)
    elseif ntens == 4 % EVZ
        % TODO: noch implementieren
        msg = 'EVZ noch nicht implementiert';
        error(msg)
        
    end
    
    sigv = sqrt( 0.5 * (...
             (sig(1)-sig(2))^2+(sig(1)-sig(3))^2+ (sig(2)-sig(3))^2 + ...
              6*(sig(4)^2+sig(5)^2+sig(6)^2) ...
                       )...
               );
    SIGV(i) = sigv;
    epsev = sqrt(0.5 * ((epse(1)-epse(2))^2+(epse(1)-epse(3))^2+ ... 
            (epse(2)-epse(3))^2+3/2*(epse(4)^2+epse(5)^2+epse(6)^2)))/(1+nu);
    epspv = sqrt(2/3*epsp'*P_tilde*epsp);
    EPSV(i) = epsev+epspv;


end