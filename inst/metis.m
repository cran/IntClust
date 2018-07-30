% function labels=metis(x,k)
%
% copyright (c) 1998-2011 by Alexander Strehl

function labels=metis(x,k,fn) 

filename = wgraph(x,[],0,fn);
labels = sgraph(k,filename,fn);
delete(filename);
