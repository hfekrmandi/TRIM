function result = iterative_ci(G,method_,eps_)
global S1
for i_agent=1:size(G.Adj,2)
    result.consenus{i_agent}.P_prior{1} =S1(:,:,i_agent);
end

for i_consensus = 2:150
    for j_agent = 1 : size(G.Adj,2)
        result.consenus{j_agent}.P_prior{i_consensus} = zeros(size(S1(:,:,1)));
        idx_neighbours = find(G.Adj(j_agent,:));
        S1_local=[];
        for j_neigh=1:length(idx_neighbours)
            S1_local(:,:,j_neigh) = result.consenus{idx_neighbours(j_neigh)}.P_prior{i_consensus-1};
        end
        result.consenus{j_agent}.P_prior{i_consensus} = (global_ci(S1_local,method_,eps_));
    end
end
