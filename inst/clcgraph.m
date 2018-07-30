% function cl = clcgraph(x,k,sfct)
%
% DESCRIPTION
%   provides cluster labels 1 to k from edge weighted value balanced
%   graph partitioning
% BUGS 
%   pmetis doesn't seem to work on PCWIN for vertex weighted graphs
%
% Copyright (c) 1998-2011 by Alexander Strehl


function cl = clcgraph(x,k,sfct,fn)

cl = cmetis(checks(simbjac(x)),sum(x,2),k,fn);
