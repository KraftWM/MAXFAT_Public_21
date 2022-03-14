function [esig,eepsp] = uniax_neuber(E,sig,epsp)
% einachsiger Neuber !!!!!  zur Ermittlung von BAUTEILFLIEßKURVEN !!!!!
%
% INPUT:
% E    -> Emodul
% sig  -> ELATISCH PLASTISCHE Spannungen
% epsp -> PLASTISCHE Dehnungen
%
% OUTPUT:
%  esig   -> pseudo elastische Spannung für esig-approach
%  eepsp  -> pseudo plastische Dehnung für eeps-approach
%
% -------------------------------------------------------------------------
esig = sqrt(sig .* (sig + E .* epsp));
eepsp = (esig-sig)./E;
end
