function [weights_ci, inf_mat, inf_vect] = calc_ci_weights_ver3(S1, local_inf_vec, method_)
    % Number of covariance matrices 
    nCovSamples = size(S1, 3);

    % Generate a random initialize weight and normalize it so 
    % it sums up to 1.
    x0 = rand(nCovSamples, 1);
    x0 = x0 ./ sum(x0);

    % Thos constraint ensures that the sun of the weights is 1
    Aeq = ones(size(x0))';
    beq = 1;

    % Weights belong to the interval [0,1]
    lb = zeros(size(x0))';
    ub = ones(size(x0))';
    A = [];
    b = [];
    nonlcon = [];

    if verLessThan('matlab', '8.4')
        options = optimset('Algorithm', 'sqp');

        if strcmp(method_, 'tr')
            x = fmincon(@cost_ci_tr, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
        elseif strcmp(method_, 'det')
            x = fmincon(@cost_ci_det, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
        end

    else
        options = optimoptions('fmincon', 'Display', 'none', 'Algorithm', 'sqp');

        if strcmp(method_, 'tr')
            x = fmincon(@cost_ci_tr, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
        elseif strcmp(method_, 'det')
            x = fmincon(@cost_ci_det, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
        end

    end

    % CI weghts
    weights_ci = x;

    % Normalize just in case
    if sum(weights_ci) > 1
        weights_ci = weights_ci ./ sum(weights_ci);
    end

    % Now that we have the weights, calculate w1*I1+...+wn*In
    inf_vect = special_dot_sum(weights_ci, local_inf_vec, 0);
    inf_mat = calc_inf_ci(x);

    % Trace cost function as the objective function
    function cost_tr = cost_ci_tr(x)
        information_matrix = zeros(size(S1(:, :, 1)));

        for i_tr = 1:length(x)
            information_matrix = information_matrix +x(i_tr, 1) * (S1(:, :, i_tr));
        end

        % Make the information matrix symetric in case numerical errors during the summation calculation
        information_matrix = 0.5 * (information_matrix + information_matrix');

        cost_tr = trace(inv(information_matrix));
    end

    % Determinant cost function
    function cost_det = cost_ci_det(x)
        information_matrix = zeros(size(S1(:, :, 1)));

        for i_det = 1:length(x)
            information_matrix = information_matrix +x(i_det, 1) * (S1(:, :, i_det));
        end

        % Make the information matrix symetric in case numerical errors during the summation calculation
        information_matrix = 0.5 * (information_matrix + information_matrix');

        cost_tr = -log(det(information_matrix));

        % cost calculation near the singularity.
        if isinf(cost_tr)
            cost_tr = log(det(inv(information_matrix)));
        end

    end

    function [information_matrix] = calc_inf_ci(x)
        information_matrix = zeros(size(S1(:, :, 1)));

        for i_det = 1:length(x)
            information_matrix = information_matrix +x(i_det, 1) * (S1(:, :, i_det));
        end

        % Make the information matrix symetric in case numerical errors during the summation calculation
        information_matrix = 0.5 * (information_matrix + information_matrix');
    end

end
