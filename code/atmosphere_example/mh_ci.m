function result = mh_ci(G)
global S1

for i_agent=1:size(G.Adj,2)
    result.consenus{i_agent}.Y_prior{1} = inv(S1(:,:,i_agent));
end

for i_consensus = 2:70
    for j_agent = 1 : size(G.Adj,2)
        result.consenus{j_agent}.Y_prior{i_consensus} = zeros(size(S1(:,:,1)));
        for k_agent = 1 :  size(G.Adj,2)
            p_jk = G.p(j_agent,k_agent);
            %             b_jk = covariance_intersection()
            updated_Y_prior = result.consenus{j_agent}.Y_prior{i_consensus} + p_jk*result.consenus{k_agent}.Y_prior{i_consensus-1};
            result.consenus{j_agent}.Y_prior{i_consensus} = updated_Y_prior;
        end
    end
end
end