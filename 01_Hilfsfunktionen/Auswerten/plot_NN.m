function ax = plot_NN(ax,Nexp,Ncalc,titleName,legendentry,...
                      handlevisible,varargin)
% Plotet N-N Diagramm
% INPUT:
% ax            - axis handle
% Nexp          - Experimentelle Ergebnisse
% Ncalc         - Rechenergebnisse 
% titlename     - Name der figure
% legendentry   - Eintrag für Legende
% handlevisible - Anzeige in Legende ('on' oder 'off')
% varargin      - plot optionen
%--------------------------------------------------------------------------

% öffne Figure
if ax == 0
    figure; grid on, hold on
    ax = gca;
    latexinterpreter;
end

% Filtern Durchläufer
NDl = 1e9;
idxCalcDl = Ncalc > NDl;
idxExpDl = Nexp > NDl;
Ncalc(idxCalcDl) = NDl;
Nexp(idxExpDl) = NDl;

% Plote N-N
plot(ax,Nexp,Ncalc,'handlevisibility',handlevisible,varargin{:});

% Plote Pfeile an die Durchläufer
% pfeil = char(8599);                                                                  
% if any(idxExpDl)
%     text(Ncalc(find(idxExpDl,1,'first')),Nexp(find(idxExpDl,1,'first')),'$\rightarrow$','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',20);
%     text(Ncalc(find(idxExpDl,1,'first')),Nexp(find(idxExpDl,1,'first')),'$\uparrow$','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',20);
% end
% logScale
set(gca,'XScale','log','YScale','log')

% namen
xlabel('$N_{exp}$'), ylabel('$N_{calc}$')
title(titleName);

% Faktor 3
% plot([1,1e10],[3,3e10],'r','LineWidth',1.3,'handlevisibility','off')
% plot([1,1e10],[0.3333,0.3333*1e10],'r','LineWidth',1.3,'handlevisibility','off')
% 
% % Faktor 10
% plot([1,1e10],[10,10e10],'b','LineWidth',1.3,'handlevisibility','off')
% plot([1,1e10],[0.1,0.1*1e10],'b','LineWidth',1.3,'handlevisibility','off')
% 
% % Mitteldiagonale
% plot([1,1e10],[1,1e10],'k','LineWidth',1.3,'handlevisibility','off')

% Achsen beschneiden
axis([1 1e7 1 1e7])

% Legende Hinzufügen
if strcmp(ax.Legend.String,'data1')
    ax.Legend.String = {legendentry};
elseif strcmp(handlevisible,'on')
    ax.Legend.String{length(ax.Legend.String)} = legendentry;
end
end