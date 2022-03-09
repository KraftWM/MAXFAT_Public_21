function latexinterpreter
% Funktion Setzt Eigenschaften der aktuellen Figure
set(gca,'FontSize',16)  
set(gca,'TickLabelInterpreter','Latex')
set(gca,'defaulttextinterpreter','latex')
set(gcf,'defaulttextinterpreter','latex')
set(legend,'Interpreter','Latex')
set(legend,'FontSize',11)
