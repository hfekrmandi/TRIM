function result = iterative_ci(G,method_)
global S1
for i_agent=1:size(G.Adj,2)
    result.consenus{i_agent}.P_prior{1} =S1(:,:,i_agent);
end

for i_consensus = 2:70
    for j_agent = 1 : size(G.Adj,2)
        result.consenus{j_agent}.P_prior{i_consensus} = zeros(size(S1(:,:,1)));
        idx_neighbours = find(G.Adj(j_agent,:));
        S1_local=[];
        for j_neigh=1:length(idx_neighbours)
            S1_local(:,:,j_neigh) = result.consenus{idx_neighbours(j_neigh)}.P_prior{i_consensus-1};
        end
        result.consenus{j_agent}.P_prior{i_consensus} = (global_ci(S1_local,method_));
    end
end


%    function cost_tr = cost_ci_tr(x)
% 
%         information_matrix = zeros(size(S1_local(:,:,1)));
%         for i_tr=1:length(x)
%             information_matrix = information_matrix +x(i_tr,1)*inv(S1(:,:,i_tr));
%         end
%         cost_tr = trace(inv(information_matrix));
%     end
%     function cost_tr = cost_ci_det(x)
%         information_matrix = zeros(size(S1_local(:,:,1)));
%         for i_det=1:length(x)
%             information_matrix = information_matrix +x(i_det,1)*inv(S1(:,:,i_det));
%         end
%         cost_tr = det(inv(information_matrix));
%     end





end