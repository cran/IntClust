% function cl = cspa(cls,k)
%
% DESCRIPTION
%  Performs CSPA for CLUSTER ENSEMBLES
%
% Copyright (c) 1998-2011 by Alexander Strehl

function cl = cspa(cls,fn,k)

disp('CLUSTER ENSEMBLES using CSPA');

if ~exist('k'),
 k = max(max(cls));
end;

clbs = clstoclbs(cls);
s = clbs' * clbs;

s = checks(s./size(cls,1));
%csvwrite('SimM.csv', s)
cl = metis(s,k,fn);

