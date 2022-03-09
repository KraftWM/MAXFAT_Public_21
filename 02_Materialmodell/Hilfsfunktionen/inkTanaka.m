function [dA,dC] = inkTanaka(ntens,A,C,n,ca,ct,dp)
% Funktion berechnet Inkrement des Nichtproportionalitätsparameters sowie
% des Tanaka Tensors
%
% INPUT:
%  ntens   - Anzahl Tensorkomponenten
%   A      - Wert NP Parameter
%   C      - Tanka Tensor (als Vektor gespeichert)
% n,nTrans - Normalentensor an FF ( und Transponierte des Normalentensors)
% ca,ct    - Modellparameter geschwindigkeit Tanaka Parameter (ca) & 
%            Tanka Tensor(ct) 
%  dp      - inkrement plastische Bogenlänge
%
% OUTPUT:
%  dA      - Inkrement NP Parameter
%  dC      - Inkrement Tanaka Tensor (als vektor)
%
% NOTATION
%
%    Tanaka Tensor
%     -> Symmetrischer Tensor 4. Stufe (Dehnungscharakter)
%     -> 21 unabhängige Komponenten (bei ntens = 6)
%     -> Speicherrichtung entlang der Diagonalen
%            Tensorkomp. in Voigt Notation                      Für Zustandsvariablen
%        ( 1111   1122   1133  2*1112  2*1113  2*1123 )       ( 1  7  12  16  19  21)
%        (        2222   2233  2*2212  2*2213  2*2223 )       (    2   8  13  17  20)
%   CT = (               3333  2*3312  2*3312  2*3323 )  =    (        3   9  14  18)
%        (                     4*1212  4*1213  4*1223 )       (            4  10  15)
%        (                             4*1313  4*1323 )       (                5  11)
%        (                                     4*2323 )       (                    6)
%
%
%__________________________________________________________________________


% -------------------------------------------------------------------------
% Unterscheide Spannungszustände
if ntens == 6 % 3D
    
    % ... C_ijkl : C_ijkl (abfangen doppelter und vierfacher Werte)
    dotC =  C(1)^2 +  C(2)^2 +  C(3)^2 + ...
            C(9)^2 + C(13)^2 + C(14)^2 + ...
           C(16)^2 + C(17)^2 + C(18)^2 + ...
           C(19)^2 + C(20)^2 + C(21)^2 + ...
       2*(  C(7)^2 +  C(8)^2 + C(12)^2 ) + ...
     0.5*( C(10)^2 + C(11)^2 + C(15)^2 ) + ...
    0.25*(  C(4)^2 +  C(5)^2 +  C(6)^2 );

    % ... n_kl * C_ijkl
    Cn = zeros(6,1);
    Cn(1) =      C( 1)*n(1) + C( 7)*n(2) + C(12)*n(3) + ...
           0.5*( C(16)*n(4) + C(19)*n(5) + C(21)*n(6) );
       
    Cn(2) =      C( 7)*n(1) + C( 2)*n(2) + C( 8)*n(3) + ...
           0.5*( C(13)*n(4) + C(17)*n(5) + C(20)*n(6) );
       
    Cn(3) =      C(12)*n(1) + C( 8)*n(2) + C( 3)*n(3) + ...
           0.5*( C( 9)*n(4) + C(14)*n(5) + C(18)*n(6) );
       
    Cn(4) =      C(16)*n(1) + C(13)*n(2) + C( 9)*n(3) + ...
           0.5*( C( 4)*n(4) + C(10)*n(5) + C(15)*n(6) );
       
    Cn(5) =      C(21)*n(1) + C(20)*n(2) + C(18)*n(3) + ...
           0.5*( C(15)*n(4) + C(11)*n(5) + C( 6)*n(6) );
       
    Cn(6) =      C(19)*n(1) + C(17)*n(2) + C(14)*n(3) + ...
           0.5*( C(10)*n(4) + C( 5)*n(5) + C(11)*n(6) );
    % ... n_kl*C_ijkl : C_ijmn*n_mn
    nCCn =         Cn(1)^2 + Cn(2)^2 + Cn(3)^2 + ...
           0.5 * ( Cn(4)^2 + Cn(5)^2 + Cn(6)^2 );
       
elseif ntens == 3 % ESZ
    
    % C_ijkl:C_ijkl
    dotC = 2 * ( ...
                   2 * ( C(1)^2 + C(2)^2 ) ...
                 + 5 * C(4)^2 ...
                 +     C(5)^2 + C(6)^2 ...
                 + 0.125 * C(3)^2 ...
                 + 4 * ( C(1)*C(4) + C(2)*C(4) ) ...
                 + C(1)*C(2) + C(5)*C(6) ...
               );
           
    % C_ijmn:n_nm
    Cn = zeros(3,1);
    Cn(1) = 2 * ( C(1) * n(1) + C(4) * n(2) ) ...
                + C(1) * n(2) + C(4) * n(1) ...
          + 0.5 * C(6) * n(3);
      
    Cn(2) = 2 * ( C(4) * n(1) + C(2) * n(2) ) ...
                + C(4) * n(2) + C(2) * n(1) ...
          + 0.5 * C(5) * n(3);
      
    Cn(3) = 2 * ( C(6) * n(1) + C(5) * n(2) ) ...
                + C(6) * n(2) + C(5) * n(1) ...
          + 0.5 * C(3) * n(3);
      
    % n_kl:C_ijkl_C_ijmn:n_mn
    nCCn = 2 * ( Cn(1)^2 + Cn(2)^2 + Cn(1)*Cn(2)) ...
           + 0.5 * Cn(3)^2;

elseif ntens == 2 % sig-tau
    
    % C_ijkl:C_ijkl
    dotC = 2.25 * C(1)^2 + 0.25 * C(2)^2 + 1.5*C(3)^2;
    
    % C_ijmn:n_nm 
    Cn = zeros(2,1);
    Cn(1) = 3 * C(1)*n(1) + C(3)*n(2);
    Cn(2) = 3 * C(3)*n(1) + C(2)*n(2);
    
    % n_kl:C_ijkl_C_ijmn:n_mn
    nCCn = 0.125 * ( 3*Cn(1)^2 + Cn(2)^2);
    
else
    msg = 'Spannungszustand für Tanaka Parameter nicht implementiert';
    error(msg)
end % Ende verzweigung Spannungszustände

% ... Target Wert np para, AT = sqrt(1 - n_kl:C_ijkl_C_ijmn:n_mn/% C_stop:C_stop ) 
AT = 0;
if dotC > 0
    AT = 1 - nCCn/dotC;
    if AT >= 0
        AT = sqrt(AT);
    else
        AT = 0;
    end
end

% ... Inkrement NP Para
dA = ca * ( AT - A ) * dp;

% ... n_ij n_kl
N = symmat2vec(n*n');
% ... ct * (N_ijkl - C_ijkl) *dp
dC = ct * (N - C)*dp;

end % Ende Funktion