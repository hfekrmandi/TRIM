% Create a k-regular Harary graph
% INPUTs: n - # nodes, k - degree of each vertex
% OUTPUTs: adj - the adjacency matrix of k-regular directed Harary graph
% AJ, Last updated: October 6, 2013
% Sample run: adj = kregularHarary(8,4)

%% function kregularHarary - main function in this file.
function adj = kregularHarary(n,k)

adj=[];

%preliminary check
if k>n-1; fprintf('a simple graph with n nodes and k>n-1 does not exist\n'); return; end
if mod(k,2)==0
    adj = hararyEven(n,k);
else
    adj = hararyEven(n,k-1);
    adj = hararyOdd(n,adj);
end
% Get only the lower triangular matrix out of adj
% this is required to prevent Matlab from creating double edges.
adj_tri = tril(adj);
bg2=biograph(adj_tri);
bg2.layoutType = 'equilibrium';
bg2.showarrows = 'off';
view(bg2);
end

%% function hararyOdd (local function)
% updates the adjacency matrix for a Harary graph for which k is odd
% returns adj - the updated adjacency matrix
function adj = hararyOdd(n,adj)
if mod(n,2) == 0
    halfDist = n/2;
else
    halfDist = (n+1)/2;
end
for i=1:n
    nextNode = i+halfDist;
    if nextNode > n
        nextNode = nextNode - n;
    end
    adj(i,nextNode)=1;
end
end

%% function harayEvn (local function)
% creates a harary graph for even degree regularity
% reutrns the adjacency matrix for the Haray Graph (n,k) where k is even
function adj = hararyEven(n,k)
adj=zeros(n,n);
halfK = k/2;
fprintf('halfK value : %d',halfK);
for i=1:n
    for j=1:halfK
        nextNode = i+j;
        if (nextNode > n)
            nextNode = nextNode - n;
        end
        prevNode = i-j;
        if (prevNode <= 0)
            prevNode = prevNode + n;
        end
        adj(i, nextNode) = 1;
        adj(i, prevNode) = 1;
    end
end

end