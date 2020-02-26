function [str]=func2strn(fh)
% function converts function handle to string containing the function name
fhstr=func2str(fh);
i1=fhstr==')'; i1i=i1.*[1:length(i1)];
ib=min(i1i(i1)); if isempty(ib); ib=0; end
i1=fhstr=='('; i1i=i1.*[1:length(i1)];
ie=max(i1i(i1)); if isempty(ie); ie=length(fhstr)+1; end
if ib>=ie; error(' Erroneous function handle'); end
str=fhstr(ib+1:ie-1);
