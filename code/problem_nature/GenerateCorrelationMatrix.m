%%
% Generates a correlation matrix with uniformly distributed correlation coefficients
% @param 	dim		Dimension of correlation matrix to generate
% @return 			Correlation matrix of dimension dim
%
function [R] = GenerateCorrelationMatrix(dim) 
	assert(dim > 0);
% 	reset(RandStream.getGlobalStream,sum(100*clock));
	
	R = ones(dim);
    L = 1;
    for i=1:dim-1
        %LInv = lu(L);
        s = SampleCorrelation(dim, i);
        r = L*s;
        
        R(i+1,1:i) = r';
        R(1:i,i+1) = r';
        L = [L zeros(i,1) ; s' sqrt(1 - s'*s)];
    end
end

%
% Calculate the distribution for a sample in dimension length(R)+1 when the
% correlation matrix R is given. 
%
function [s] = SampleCorrelation(N, i)
    m = (N - i - 1)/2;
    
    % random sample on n-sphere
    u=randn(i, 1);
    u=u/norm(u);
    
    % determine distance from 0
    % x = betarnd(i/2, m+1);
    x = BetaRandFast(i/2, m+1);
    x = sqrt(x);
    
    % gives point in sphere
    s = u*x;
    
    % give back ellipsoidal point
    %r = L * s;
end

function ret = BetaRandFast(a,b)
    x = randg(a, 1);
    ret = x ./ (x + randg(b,1));
end