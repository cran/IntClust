% function labels=cmetis(e,w,k)
%
% copyright (c) 1998-2011 by Alexander Strehl


function labels=cmetis(e,w,k,fn) 

filename = wgraph(e,w,1,fn);
labels = sgraph(k,filename,fn);
delete(filename);
