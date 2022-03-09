function [G, P] = esedableitung(SIG,EPS,...
                                  REFSIG,REFEPS)

% =====================================================================
% P
P = diag(SIG - REFSIG);

% =====================================================================
% G
G = zeros(3,3);

end % Ende Funktion