function OPTIONS=constopt(parain)
%constopt Default parameters for use with CONSTR.
%
%   OPTIONS = constopt returns the default options.
%
%   OPTIONS = constopt(SOME_OPTIONS) takes the non-empty vector 
%     SOME_OPTIONS and replaces the zero or missing values of 
%     SOME_OPTIONS with the default options.
%
%   The parameters are:
%   OPTIONS(1)-Display parameter (Default:0). 1 displays some results.
%   OPTIONS(2)-Termination tolerance for X.(Default: 1e-4).
%   OPTIONS(3)-Termination tolerance on F.(Default: 1e-4).
%   OPTIONS(4)-Termination criterion on constraint violation.(Default: 1e-4)
%   OPTIONS(5)-Algorithm: Strategy:  Not always used.
%   OPTIONS(6)-Algorithm: Optimizer: Not always used. 
%   OPTIONS(7)-Algorithm: Line Search Algorithm. (Default 0)
%   OPTIONS(8)-Function value. (Lambda in goal attainment. )
%   OPTIONS(9)-Set to 1 if you want to check user-supplied gradients
%   OPTIONS(10)-Number of Function and Constraint Evaluations.
%   OPTIONS(11)-Number of Function Gradient Evaluations.
%   OPTIONS(12)-Number of Constraint Evaluations.
%   OPTIONS(13)-Number of equality constraints. 
%   OPTIONS(14)-Maximum number of function evaluations. 
%               (Default is 100*number of variables)
%   OPTIONS(15)-Used in goal attainment for special objectives. 
%   OPTIONS(16)-Minimum change in variables for finite difference gradients.
%   OPTIONS(17)-Maximum change in variables for finite difference gradients.
%   OPTIONS(18)-Step length. (Default 1 or less). 
%   OPTIONS(19)-Flag for succesful termination 1 if successful %GvW 2-5-2007
%   OPTIONS(20)-Merit function factor, amplifies function value w.r.t. constraints %GvW
%   OPTIONS(21)-CPU time %GvW
%   OPTIONS(22)-Max(g), tolerance by which the solutions satisfies the constraints
%   OPTIONS(23)-Maximum number of consecutive small steps
%   OPTIONS(24)-Maximum number of consecutive failures of constrained & function improvement
%   OPTIONS(25)-Maximum number of consecutive failures of constrained satisfaction


%   Copyright 1990-2003 The MathWorks, Inc. 
%   $Revision: 1.1.6.1 $  $Date: 2003/12/31 19:52:53 $

if nargin<1; 
   parain = [];
end
sizep=length(parain);
OPTIONS=zeros(1,25); % GvW 2-5-2007 Modified from zeros(1,18)
OPTIONS(1:sizep)=parain(1:sizep);
default_options=[0,1e-4,1e-4,1e-4,0,0,0,0,0,0,0,0,0,0,0,1e-8,0.1,0,0,0,0,0,0,0,0];
OPTIONS=OPTIONS+(OPTIONS==0).*default_options;
