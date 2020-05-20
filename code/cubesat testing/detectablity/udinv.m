% UDINV : inverse of unit upper triangular u and diagonal d
%
%         function [udi,r]=udinv(ud,tolp,toln)
%
%	  udi : inv(ud), d on diagonal
%         r   : rank(di)
%         tolp: tolerance for positiveness
%               default 1e-12, see also toln
%         toln: error generation if
%	        diagonal part < toln
%               default: tolp
%               if toln < 0 all the diagonal
%	        parts < max(0,tolp) are set to zero
%
  function [udi,r]=udinv(ud,tolp,toln)

  if nargin > 3; error('  one, two or three input arguments required'); end;
  if nargin==1; tolp=1e-12; toln=tolp;
  elseif nargin==2; toln=tolp; end; tolp=max(0,tolp); 
  if (toln>tolp); error('  toln > tolp'); end;

  [n,m]=size(ud);
  if n~=m; error(' ud must be square'); end;
  if n==0; error('  Compatible but empty inputs'); end;

  r=n;
  for i=1:n
    if ud(i,i)<toln; error('  toln violated');
    elseif ud(i,i)<=tolp; udi(i,i)=0; r=r-1;
    else; udi(i,i)=1/ud(i,i); end
  end;

  for j=2:n
    jm1=j-1;
    for k=1:jm1
      sum=-ud(k,j);
      for i=k+1:jm1
        sum=sum-udi(k,i)*ud(i,j);
      end
      udi(k,j)=sum;
    end
  end
