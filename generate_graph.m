function G = generate_graph(Adj)
    % This function accepts an adjecancy matrix where degree of each
    % node is equal to 1 + number of its neighbours. That is, all agents
    % are connected to themselves as well.

    %  number of nodes ( = |v|)
    nv = size(Adj, 2);

    % Assign the graph adjecancy matrix
    G.Adj = Adj;

    % Calculate inclusive node degrees
    for i = 1:nv
        G.d(i) = sum(G.Adj(i, :)) + 1;
    end

    % Calculate weights for MHMC distributed averaging
    % This is slightly different from the formula used in the paper
    % http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.161.3893&rep=rep1&type=pdf
    for i = 1:nv

        for j = 1:nv

            if G.Adj(i, j) ~= 0

                if i ~= j
                    G.p(i, j) = min(1 / G.d(i), 1 / G.d(j));
                end

                if i == j
                end

            else
                G.p(i, j) = 0;
            end

        end

    end

    for i = 1:nv
        G.p(i, i) = 1 - sum(G.p(i, :));
    end

end
