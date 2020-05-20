function [tu]=pwaver(tu,t)
% PWAVER: Piecewise average u from tu table based on t
%
%         [tu]=pwaver(tu,t)
%
%         Input:
%         tu   :  [t u] table
%         t    :  Time axis for piecewise averaging u. t must be
%                 inside the range tu(:,1) and monotonically increasing
%
%         Output:
%         tu   :  Piecewise averaged [t u] table
%
% GvW 21-01-2010

% Check inputs
if nargin~=2; error(' Two inputs required'); end

% Sizes
[n,m]=size(t);
if m~=1;
  if n~=1; error(' t must be a vector'); else t=t'; n=m; end
end
[n1,m1]=size(tu);

% Obtain time range
tmax=max(tu(:,1)); tmin=min(tu(:,1));

% Interpolate tu with a 10 times finer grid than t
dt=(t(end)-t(1))/(n-1); tn=[t(1):dt/10:t(end)]';
u=interp1(tu(:,1),tu(:,2:end),tn); tu=[tn u];

% Compute averages u from new tu
u=zeros(n,m1-1);
for j=1:n-1
  if t(j)<tmin || t(j+1)>tmax
    error(' t out of range tu(:,1)')
  elseif t(j+1)<=t(j)
    error(' t not monotonically increasing')
  end
  % Get u values in between t(j) and t(j+1)
  ih=tu(:,1)>=t(j) & tu(:,1)<=t(j+1);
  if max(ih)==0;
    u(j,:)=interp1(tu(:,1),tu(:,2:end),t(j)); % Take interpolated values
  else    
    u(j,:)=mean(tu(ih,2:end)); % Take mean values
  end
end

% Final u & generation new tu
u(n,:)=zeros(1,m1-1); tu=[t u];