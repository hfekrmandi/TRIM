function [y]=iscomp(s)

c=[{'<'} {'<='} {'>'} {'>='} {'=='}]; 
if ~isstr(s) | ~strcmp(s,c); y=0; else y=1; end