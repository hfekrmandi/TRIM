function [perf_index_prob,perf_index] = atmosphere_diagram_shell()
clear all
global opt_dist
profile on
close all
% dbclear all
dbstop if error
dbstop if warning
[A,x0,B,C] = create_sys_atmosphere();



flag_converged = 0;
global fail_prob reg_deg
range_prob = [ 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
% range_prob = [ 0.4 0.6 0.8  1];
% range_prob = [ 0 0.2  0.4 0.6 0.8 1];
range_reg = [  4];

h_res = figure;
% range_prob = [ 1];
for j_reg=1:length(range_reg)
    reg_deg = range_reg(j_reg);
    for i=1:length(range_prob)
        
        fail_prob =  range_prob(i);
        problem_def_diag1(A,x0(1:size(A,1)));
        opt_dist.reg_degree = reg_deg;
        tic
        for i_step = 1:20
            i_step
            opt_dist.i_step = i_step;
            %     flag_converged = 0;
            sim_system();
            pred();
            consenus_ver2();
            %         update_();
            [perf_index(i_step,:),error_results{i_step},opt_dist_result{i_step},bcs_perf(i_step,:)] = post_process_diag_ver2();
            %         if opt_dist.FLAGS.debug
            %             diagnostic_tests(error_results{i_step},opt_dist_result{i_step})
            %         end
        end
        profile viewer
        time_(i) = toc;
        snapshot_perf{i} = perf_index;
        snapshot_perf_bcs{i} = bcs_perf;
        perf_index_prob(j_reg,i,:) = mean(bcs_perf);
        %     save(['file_0',num2str(i)],'bcs_perf','error_results', 'perf_index', 'range_prob', 'x0')
        %     save(['results_',num2str(problem_def_diag1)]);
        %     save(['atmosh_reduced_order',num2str(i)],'bcs_perf','error_results', 'perf_index', 'range_prob', 'x0')
        %     save(['atmosh_shell',num2str(i),'_',num2str(j_reg)],'bcs_perf','error_results', 'perf_index', 'range_prob', 'x0')
        
    end
    disp('aa')
end
fclose('all')
% save('example_agents_CI_Comp_70_steps_4')
close all
figure
range_prob =[ 0.4 0.6 0.8  1];

plot(range_prob,100*perf_index_prob(:,1)','LineWidth',2); hold on
plot(range_prob,100*perf_index_prob(:,2)','LineWidth',2); hold on
legend('our method','CI')
xlabel('Link Fail Prob.')
ylabel('Est.  Performance')
grid on
grid minor

figure

plot(perf_index(:,1)); hold on
plot(perf_index(:,2)); hold on;
end

