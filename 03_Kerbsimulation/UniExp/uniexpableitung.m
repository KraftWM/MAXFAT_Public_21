function [G, P] = uniexpableitung(SIG,EPS,...
                                  REFSIG,REFEPS,Cq)



% =====================================================================
% P
P = (1+Cq) * diag(SIG - REFSIG);

% =====================================================================
% G
G = (1-Cq) * diag(EPS - REFEPS);

end % Ende Funktion