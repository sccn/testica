function [retval] = relvar(d1, d2);

% RELVAR returns the relative residual variance between measured and simulated data
% 
% rv = relvar(measured, simulated)

% Copyright (C) 1999, Robert Oostenveld
%
% $Log: relvar.m,v $
% Revision 1.2  2003/03/11 14:45:37  roberto
% updated help and copyrights
%

retval = sum((d1-d2).^2) ./ sum(d1.^2);

