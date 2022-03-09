function [P, P_line, P_hat, A, P_check] = set_maps(ntens,ndi)
% Funktion setzt Abbildungsmatrizen je nach Spannungszustand 
%
% INPUT:
%    ntens -> Anzahl Tensorkomponenten
%    ndi   -> Anzahl Diagonalelemente
% OUTPUT:
%  Allgemein
%    P     -> Abbildung tensor auf deviator
%  P_line  -> Hilfsmatrix Skalaprodukte zwischen Spannungen (dev Span in
%                                                            ESZ)
%  speziell ESZ & sig-tau
%   P_hat  -> Skalarprodukte Spannungen (min eine span kein dev)
%             und Dehnungen aus Spannunge
%     A    -> Spannungen aus Dehnungen 
%  P_check -> Skalarprodukte dev Spannungen und Dehnungen (ESZ)
%__________________________________________________________________________

% -------------------------------------------------------------------------
%                    3D
if ntens == 6 % 3D
    
    P = 1/3 .* [ 2,-1,-1, 0, 0, 0;...
                -1, 2,-1, 0, 0, 0;...
                -1,-1, 2, 0, 0, 0;...
                 0, 0, 0, 3, 0, 0;...
                 0, 0, 0, 0, 3, 0;...
                 0, 0, 0, 0, 0, 3];    
             
    P_line = diag([1,1,1,2,2,2]);   
    
    if nargout == 3
        P_hat = diag([1,1,1,0.5,0.5,0.5]); 
    elseif nargout == 4
        P_hat = diag([1,1,1,0.5,0.5,0.5]);
        A = diag([1 1 1 1 1 1]);
    elseif nargout == 5
        P_hat = diag([1,1,1,0.5,0.5,0.5]);
        A = diag([1 1 1 1 1 1]);
        P_check = diag([1 1 1 1 1 1]);
    end
% -------------------------------------------------------------------------
%                    ESZ
elseif ntens == 3 && ndi == 2 % ESZ
    
    P = 1/3 .* [ 2,-1, 0;...
                -1, 2, 0;...
                 0, 0, 3];   

    P_line = [2,1,0;...
              1,2,0;...
              0,0,2];  

    if nargout == 3
        
        P_hat = diag([1,1,2]);

        
    elseif nargout == 4 
        
        P_hat = diag([1,1,2]);

        A = diag([1,1,0.5]);

        
    elseif nargout == 5
        
        P_hat = diag([1,1,2]);

        A = diag([1,1,0.5]);

        P_check = [2,1,0;...          
                   1,2,0;...
                   0,0,1]; 

    end
% -------------------------------------------------------------------------
%                    sigma - tau
% !!! Hier alle Abbildungen als Vektor da alles Diagonalmatricen sind
elseif ntens == 2 && ndi == 1
    P = 1/3 .* [ 2;...
                 3];   

    P_line = 1/2 .* [3;...
                     4];    

    if nargout == 3
        
        P_hat = [1;2];

        
    elseif nargout == 4 
        
        P_hat = [1;2];

        A = [1;0.5];

        
    elseif nargout == 5
        
        P_hat = [1;2];

        A = [1;0.5];

        P_check = 0.5*[3;2]; 

    end
% -------------------------------------------------------------------------
%                    1D      
elseif ntens == 1
    
    P = 1;
    
    P_line = 1;
    
    P_hat = 1;
    
    % Damit Compiler im Coder nicht rumzickt

    A = 1;

    P_check = 1;
    
else
    P = 1;
    
    P_line = 1;
    
    P_hat = 1;

    A = 1;

    P_check = 1;
end