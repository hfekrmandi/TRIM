function [xo,yo] = stair(x,y,str)
%STAIR	Stairstep graph (bar graph without internal lines).
%       The difference with STAIRS is that unevenly spaced
%      	time intervals are allowed and the final time may
%       therefore be included in x (if not the final time
%       interval is assumed equal to the previous one) 
%
%	Stairstep plots are useful for drawing time history plots of
%	digital sampled-data systems.
%	STAIR(X,Y) draws a stairstep graph of the columns in matrix Y at
%	the locations specified in X.  The X-values must be in
%	ascending order but not necessarily evenly spaced.
%	[XX,YY] = STAIR(X,Y) does not draw a graph, but returns vectors
%	X and Y such that PLOT(XX,YY) is the stairstep graph.
%
%	See also BAR, HIST, STAIRS

%	L. Shure, 12-22-88.
%	Revised GvW 9-5-2001.
%	Copyright (c) 1984-94 by The MathWorks, Inc.

[nx,n] = size(x); 
if isstr(x)
	error('Input arguments must be numeric.')
end
if min(n,nx) > 1
	error('x must be a vector')
end
if nargin==3;
  if ~isstr(str);
    error(' 3rd input must be a plot string')
  end
end

% Transpose x if row vector
if (nx == 1 & n > 1), x=x'; nx=n; end 

[ny,n] = size(y);

% Transpose y if necessary
if ny~=nx & ny~=nx-1
  if n==nx | n==nx-1; y=y'; ny=n;
  else; error('x,y have incompatible dimensions'); end
end

% Create last element x if necessary
if ny==nx;  x=[x;  2*x(end)-x(end-1)]; end

% Double data except for first and last element of x
xo=x(1,:); yo = [];
for i=1:ny
  h = [ y(i,:); y(i,:)]; yo = [yo; h];
  h = [ x(i+1,:); x(i+1,:)]; xo = [xo; h];
end
xo(end,:) = [];

% Plot or output
if (nargout == 0); 
  if nargin<3
    plot(xo,yo); clear xo yo;
  else
    plot(xo,yo,str); clear xo yo;
  end;
end