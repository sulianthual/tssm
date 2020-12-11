function [fout]=ensospe_ncdfgetvar(filenc,varnc)
% by Sulian Thual
fout=ncread(filenc,varnc);
