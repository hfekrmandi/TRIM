function [m,ind]=stackomp(m1,opstr,m2,wo)
% STACKOMP : Function compares the stacks of the inputs using the operator
%            string. It returns the elements and their indices if true
%
%            Input:
%            m1, m2: two matrices of equal stack dimension
%                    m2 possibly scalar.
%            opstr : operator string 
%            wo    : if nonzero unequal size warning not displayed
%                    (optional)
%         
%            Output
%            m     : column vector with true entries m1
%            ind   : stack indices of true entries m1
%
%  Examples :
%  A=[1 3; 3 4]; B=3; [m,ind]=stackomp(A,'==',B)
%  gives: m=[3;3]; ind=[2;3].
%
%  A=[1 2; 3 5]; B=[1 3; 2 4]; [m,ind]=stackomp(A,'>',B)
%  gives: m=[3;5], ind=[2;4].
%
%  To set all matrix elements of A equal to 3 to 1:
%  [m,ind]=stackomp(A,'==',3); As=A(:); As(ind)=1; A(:)=As;
%  Alternatively: A=A-2*(A==3)
%
% GvW 8-4-2011

% Check inputs
if nargin<3; error(' 3 or 4 inputs required');
elseif nargin<4; wo=0; end
  
if ~isstr(opstr) || ~iscomp(opstr)
  error(' 2nd argument must be <, >, <=, >=, or ==');
end

s1=size(m1); s2=size(m2);

if max(s2)~=1
  if length(m1)~=length(m2)
    error(' 1st and 3rd input must have equal stack length');
  elseif (length(s1)~=length(s2) || max(s1-s2)~=0) && ~wo
    warning(' 1st and 3rd input are unequal in size');
  end
end

% Compare
m1=m1(:); m2=m2(:); im=[1:length(m1)]';
exstr=['ind=im.*(m1' opstr 'm2);'];
eval(exstr); ind=ind(ind>0); m=m1(ind);