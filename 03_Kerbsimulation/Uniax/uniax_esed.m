% =========================================================================
%            Hilfsfunktionen für einachsige Kerbnäherung
% =========================================================================
function [esig,eepsp] = uniax_esed(E,nstrich,sig,epsp)
% einachsige esed !!!!!  zur Ermittlung von BAUTEILFLIEßKURVEN !!!!!
%
% INPUT:
% E        -> Emodul
% nstrich  -> zyklischer Parameter aus ramberg osgood
% sig      -> ELATISCH PLASTISCHE Spannungen
% epsp     -> PLASTISCHE Dehnungen
%
% OUTPUT:
%  esig   -> pseudo elastische Spannung für esig-approach
%  eepsp  -> pseudo plastische Dehnung für eeps-approach
%
% -------------------------------------------------------------------------
fak = 1 - nstrich/(nstrich + 1);
esig = sqrt( sig.^2 + 2.* E .* fak .* sig .* epsp );
eepsp = (esig-sig)./E;
end