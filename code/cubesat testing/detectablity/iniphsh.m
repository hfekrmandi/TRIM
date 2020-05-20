% INIPHSH : Function generates initial values ph, sh
%
%           [ph,p12,p2]=iniphsh(ic,nx,nc)
%
%           Input:
%           ic    : initialization type
%           ic=1: 0.9*eye(nx)+0.1*ones(nx)
%           ic=2: random
%           ic=3: eye(nx)
%     	    nx :  dimension system state
%           nc :  dimension compensator state
%
%           Output :
%           ph=p12*pinv(p2)*p12'
%           p12 nx x nc
%           p2>=0 nc x nc
%
%           L.G. Van Willigenburg, W.L. De Koning, 27-06-96.
%
  function [ph,p12,p2]=iniphsh(ic,nx,nc)

  if nc<1; ph=zeros(nx); p12=[]; p2=[];
  else
    if ic<1.5;
      h1=0.1;
      p12=(1-h1)*eye(nx,nc)+h1*ones(nx,nc);
      p2=(1-h1)*eye(nc)+h1*ones(nc,nc);
    elseif ic<2.5;
      p12=rand(nx,nc); p2=psdr(nc); 
    else
      p12=eye(nx,nc); p2=eye(nc);
    end
    ph=p12*pinv(p2)*p12';
  end
  
% Small initial values
%  epsi=1e-6; ph=epsi*ph; p12=epsi*p12; p2=epsi*p2;
