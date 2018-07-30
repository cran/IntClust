% function labels=hmetis(x,k,w)
%
% copyright (c) 1998-2011 by Alexander Strehl

function labels=hmetis(x,k,fn,w) 

if ~exist('w'),
  filename = wgraph(x,[],2,fn);
else
  filename = wgraph(x,w,3,fn);
end; 
labels = sgraph(k,filename,fn); 
delete(filename);
