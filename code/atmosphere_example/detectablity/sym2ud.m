% SYM2UD: UDU factorization P=UDU', P>=0
%
%         See Bierman 1977: Factorization Methods for
%         discrete sequential estimation
%
%         [ud,r]=sym2ud(p,tolp,toln);
%
%         p     : square symmetric matrix >= 0 (see toln)
%         tolp  : tolerance for positiveness
%                 default 1e-12, see also toln
%         toln  : error generation if
%     	          to be square rooted part < toln
%                 default: tolp
%                 if toln < 0 all the to be square rooted
%                 parts < max(0,tolp) are set to zero
%         ud    : factorization matrix, u unit upper triangular
%                 d diagonal, both stored in ud
%         r     : rank of u and d
%
function [ud,r]=sym2ud(p,tolp,toln);

if nargin > 3; error('  one, two or three input arguments required'); end;
if nargin==1; tolp=1e-12; toln=tolp;
elseif nargin==2; toln=tolp; end; tolp=max(0,tolp); 
if (toln>tolp); error('  toln > tolp'); end;
%toln,tolp
[n,m]=size(p);
if n~=m; error('  p must be square'); end;
if n==0; error('  Compatible but empty inputs'); end;

ud=zeros(n,n); r=n;

% copy upper triangular part p to ud
for j=1:n; for i=1:j; ud(i,j)=p(i,j); end; end;

for j=n:-1:2
   if ud(j,j)<toln; error('  toln or symmetry violated');
   elseif ud(j,j)<=tolp; ud(j,j)=0; a=0; r=r-1;
   else; a=1./ud(j,j); end
   for k=1:j-1;
      b=ud(k,j); ud(k,j)=a*b; sum=0;
      for i=1:k
         ud(i,k)=ud(i,k)-b*ud(i,j);
      end
   end
end
if ud(1,1)<toln; error('  toln or symmetry violated');
elseif ud(1,1)<=tolp; ud(1,1)=0; r=r-1;
end
