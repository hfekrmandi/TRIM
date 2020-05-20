function [x,OPTIONS,lambda,HESS]=constr(FUN,x,OPTIONS,VLB,VUB,GRADFUN,varargin)
%CONSTR Finds the constrained minimum of a function of several variables.
%   
%   X=CONSTR(FUN,X0) starts at X0 and finds a constrained minimum to 
%   the function which is described in FUN (function handle/string with function name).
%   The function FUN should return two arguments: a scalar value of the 
%   function to be minimized, F, and a matrix of constraints, G: 
%   [F,G]=FUN(X).
%
%   X=CONSTR(FUN,X,OPTIONS) allows a vector of optional parameters to 
%   be defined. For more information type HELP constopt.
%
%   Equality constraints must appear first in G(:) followed by
%   the inequality constraints. Options(13) specifies the number of
%   equality constraints (default = 0). Then the constraints are as
%   follows: G(1:options(13)) = 0; G(options(13)+1:end) <= 0
%
%   X=CONSTR(FUN,X,OPTIONS,VLB,VUB) defines a set of lower and upper
%   bounds on the design variables, X, so that the solution is always in 
%   the range VLB <= X <= VUB. 
%   
%   X=CONSTR(FUN,X,OPTIONS,VLB,VUB,GRADFUN) allows a function 
%   GRADFUN to be entered which returns the partial derivatives of the 
%   function and the constraints at X:  [gf,GC] = GRADFUN(X).
%   Use OPTIONS(9)=1 to check analytic gradients in GRADFUN against
%   numeric gradients during the first iteration.
%
%   X=CONSTR(FUN,X,OPTIONS,VLB,VUB,GRADFUN,P1,P2,...) passes the 
%   problem-dependent parameters P1,P2,... directly to the functions FUN 
%   and GRADFUN: FUN(X,P1,P2,...) and GRADFUN(X,P1,P2,...).  Pass
%   empty matrices for OPTIONS, VLB, VUB, and 'GRADFUN' to use the 
%   default values.
%
%   [X,OPTIONS]=CONSTR(FUN,X0,...) returns the parameters used in the 
%   optimization method.  For example, OPTIONS(10) contains the number 
%   of function evaluations used.
%
%   [X,OPTIONS,LAMBDA]=CONSTR(FUN,X0,...) returns the Lagrange multipliers
%   at the solution X in the vector LAMBDA.
%
%   [X,OPTIONS,LAMBDA,HESS]=CONSTR(FUN,X0,...) returns the quasi-Newton
%   approximation to the Hessian matrix at the solution X.
%

%   Copyright 1990-2002 The MathWorks, Inc. 
%   $Revision: 1.51 $  $Date: 2002/03/12 20:36:09 $
%   Revisions by GvW 22-02-2010

if nargin < 2, error('constr requires two input arguments'); end
if nargin < 3, OPTIONS=[]; end
if nargin < 4, VLB=[]; end
if nargin < 5, VUB=[]; end
if nargin < 6, GRADFUN=[]; end

caller='constr';
lenVarIn = length(varargin);

if isa(FUN,'function_handle'); % function handle: convert to string
  FUNS=func2strn(FUN);
else % copy function string
  FUNS=FUN;
end
% Convert to inline function as needed.
gradflag = 0;
[funfcn,msg] = prefcnchk(FUNS,caller,lenVarIn,gradflag);
if ~isempty(msg)
  error(msg);
end

if ~isempty(GRADFUN)
  if isa(GRADFUN,'function_handle'); % function handle: convert to string
    GRADFUNS=func2strn(GRADFUN);
  else % copy function string
    GRADFUNS=GRADFUN;
  end
  gradflag = 1;
  [gradfcn,msg] = prefcnchk(GRADFUNS,caller,lenVarIn,gradflag);
  if ~isempty(msg)
    error(msg);
  end
else
  gradfcn = [];
end

if isa(FUN,'function_handle'); funfcn{1}=FUN; end
if isa(GRADFUN,'function_handle'); gradfcn{1}=GRADFUN; end
[x,OPTIONS,lambda,HESS]=nlconstold(funfcn,x,OPTIONS,VLB,VUB,gradfcn,varargin{:});

