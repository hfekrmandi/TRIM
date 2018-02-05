
test_idx = '3'
switch test_idx
    case '1'
        figure;imagesc(pinv(opt_dist.result.est{i_consensus}.Y_cen ))
        p1 = pinv(opt_dist.result.est{i_consensus}.Y_cen );
        p = p1;
        p = 0.5*(p+p');
        issymmetric(p)
        rcond(p)
        rcond(p1)
        for bigUncertainty = 0.05:0.01:1
            bigUncertainty = 0.07
            P_init_ = bigUncertainty*diag(~idx_ith_agent) + diag(0.05.*idx_ith_agent);
            rcond(P_init_);
            disp(['big_incertainty = ',num2str(bigUncertainty),'---rcond(P_init_) = ',num2str(rcond(P_init_))])
        end
        close all
        [A,x0,B,C] = create_sys_atmosphere();
        h_predict = figure;
        P_test = P_current;
        
        Q = 0.005*eye(size(P_test));
        for i=1:50
            P_test = opt_dist.A*P_test*opt_dist.A' + Q;
            rcond(P_test)
            disp(['big_incertainty = ',num2str(i),'---rcond(P_test) = ',num2str(rcond(P_test))])
            figure(h_predict);imagesc(P_test)
            pause(0.1)
        end
        
        NUM_PRIMAL          = 200;
        NUM_STEP = 10;
        [Ar, Br, Cr, VF, UF] = rpod(A, B, C, NUM_NON, NUM_PRIMAL, NUM_STEP)
        
        
        signal = opt_dist.A*x_current + opt_dist.B*( opt_dist.source.Q');
        noise = delta_noise
        
        signal_to_noise = signal./noise
        mean(signal_to_noise)
        
        rcond(opt_dist.result.pred.P_cen)
        rcond(opt_dist.result.pred.P_bar(:,:,i_agent))
    case '2'
        clear all
        close all
        clc
        global opt_dist
        [A,x0,B,C] = create_sys_atmosphere();
        problem_def_diag1(A,x0);
        rank(obsv(A,C))
        rank(obsv(A,C(1:7,:)))
        x_prior = x0;
        x_gt = x0;
        q_var = 0.001;
        r_var = 0.001;
        Q = q_var*eye(size(A));
        Q_input = q_var*eye(size(B,2));
        R = r_var*eye(size(C,1));
        p_prior = Q;
        C = C(1:3,:);
        h_predict = figure;
        h_estimate = figure;
        for i=1:400
            delta_noise =  randn(size(x_prior))*sqrt(q_var);
            x_predict = A*x_prior +  B*( opt_dist.source.Q');
            %       P_predict = A*p_prior*A' + Q;
            P_predict = A*p_prior*A' + B*Q_input*B';
            rcond_predcit(i) = rcond(P_predict);
            disp(['step = ',num2str(i),'---rcond(P_predict) = ',num2str(rcond(P_predict))])
            
            %       x_gt = A*x_gt +  B*( opt_dist.source.Q')+ delta_noise;
            x_gt = A*x_gt +  B*( opt_dist.source.Q' + sqrt(q_var) );
            
            
            %       R = 0.1.*(C*x_gt);
            R = r_var*eye(size(C,1));
            
            z = C*x_gt + randn(size(C*x_gt))*sqrt(r_var);
            %       z = C*x_gt + randn(size(C*x_gt)).*R;
            h = C*x_predict;
            H = C;
            
            %       K = P_predict*H'*inv(H*P_predict*H' + diag(R));
            K = P_predict*H'*inv(H*P_predict*H' + (R));
            
            x_est  = x_predict + K*(z-h);
            
            p_est = (eye(size(P_predict)) - K*H)*P_predict;
            p_est = 0.5.*(p_est+p_est');
            rcond_est(i) = rcond(p_est);
            p_error(i) = norm(p_est - P_predict)
            disp(['step = ',num2str(i),'---rcond(p_est) = ',num2str(rcond(p_est))])
            disp(['step = ',num2str(i),'---norm(p_est - p_prior) = ',num2str(norm(p_est - p_prior))])
            %       figure(h_predict)
            %       subplot(121); imagesc(P_predict); hold on;
            %       subplot(122); imagesc(p_est); hold on;
            %       pause(0.1)
            hold off
            
            x_prior = x_est;
            p_prior = p_est;
            
            
        end
        figure
        subplot(211)
        plot([100:400],rcond_est(100:end)); hold on
        plot([100:400],rcond_predcit(100:end)); hold on
        legend('reciprocal condition number of P_{est}')
        subplot(212)
        plot(p_error)
        
        %
        figure; plot(100.*perf_index(:,1)); hold on
        plot(100.*perf_index(:,2)); hold on
        legend('our','pure_{CI}')
        title ('perdormance comparison probablity of failure = 0.2')
        grid on
        grid minor
    case '3'
        close all
        e_bar_=[];
        e_bar_vs_cen=[];
        e_P_bar_cen=[];
        e_bar_CI=[];
        e_bar_CI_vs_cen=[];
        e_P_bar_CI_cen=[];
        e_bar_BC_dist = [];
        e_bar_CI_BC_dist = [];
        con_perc_ci = [];
        con_perc_bar = [];
        con_perc_cen = [];
        e_P_bar_cen_2 = [];
        e_P_bar_CI_cen_2 = [];
        e_P_bar_root =[]
        e_P_CI_root =[];
                range_1 = [10:50] ;

        for i=1:numel(error_results)
            e_bar_=[e_bar_;error_results{i}.e_bar_(range_1,:)];
            e_bar_vs_cen=[e_bar_vs_cen,error_results{i}.e_bar_vs_cen(range_1,:)];
            e_P_bar_cen=[e_P_bar_cen;error_results{i}.e_P_bar_cen(range_1,:)];
            e_P_bar_cen_2=[e_P_bar_cen_2;error_results{i}.e_P_bar_cen_2(range_1,:)];
            e_P_bar_root =[e_P_bar_root;error_results{i}.e_P_bar_root(range_1,:)]
            e_P_CI_root =[e_P_CI_root;error_results{i}.e_P_CI_root(range_1,:)]
            
            
            e_bar_CI=[e_bar_CI;error_results{i}.e_bar_CI(range_1,:)];
            e_bar_CI_vs_cen=[e_bar_CI_vs_cen;error_results{i}.e_bar_CI_vs_cen(range_1,:)];
            e_P_bar_CI_cen=[e_P_bar_CI_cen;error_results{i}.e_P_bar_CI_cen(range_1,:)];
            e_P_bar_CI_cen_2=[e_P_bar_CI_cen_2;error_results{i}.e_P_bar_CI_cen_2(range_1,:)];
            
            e_bar_BC_dist = [e_bar_BC_dist;error_results{i}.e_bar_BC_dist(range_1,:)];
            e_bar_CI_BC_dist = [e_bar_CI_BC_dist;error_results{i}.e_bar_CI_BC_dist(range_1,:)];
            
            con_perc_ci = [con_perc_ci; error_results{i}.con_perc_ci(range_1,:)];
            con_perc_bar = [con_perc_bar; error_results{i}.con_perc_bar(range_1,:)];
            con_perc_cen = [con_perc_cen; error_results{i}.con_perc_cen(range_1)];
            
        end
        h_err1 = figure
        h_err2 = figure
        h_err3 = figure
        range_ = [1:56] ;
        for i_agent = 1:9
            figure(h_err1)
            subplot(321)
            plot(e_bar_(:,i_agent)); hold on
            title('e_{bar} = ||x_{bar} - x_{gt}||')
            grid on; grid minor;
            
            subplot(322)
            plot(e_bar_vs_cen(:,i_agent)); hold on
            title('e_{bar vs. cen}=||x_{bar} - x_{cen}||')
            grid on; grid minor;
            
            subplot(323)
            plot(e_P_bar_cen(:,i_agent)); hold on
            title('e_{P_{bar} vs. P_{cen}}= det(P_{cen})/det(P_{bar})')
            grid on; grid minor;
            
            subplot(324)
            plot(e_bar_CI(:,i_agent)); hold on
            title('e_{CI vs. gt}=||x_{CI} - x_{gt}||')
            grid on; grid minor;
            
            subplot(325)
            plot(e_bar_CI_vs_cen(:,i_agent)); hold on
            title('e_{CI vs. cen}=||x_{CI} - x_{cen}||')
            grid on; grid minor;
            
            subplot(326)
            plot(e_P_bar_CI_cen(:,i_agent)); hold on
            title('e_{P_{CI} vs. P_{cen}}= det(P_{cen})/det(P_{CI})')
            grid on; grid minor;
            
            figure(h_err2)
            subplot(221)
            plot(e_bar_BC_dist(:,i_agent)); hold on
            title('BCdist_{bar}')
            grid on; grid minor;
            subplot(222)
            plot(e_bar_CI_BC_dist(:,i_agent)); hold on
            title('BCdist_{CI}')
            grid on; grid minor;
            subplot(2,2,[3 4])
            plot(e_bar_BC_dist(:,i_agent)); hold on
            plot(e_bar_CI_BC_dist(:,i_agent)); hold on
            grid on; grid minor;
            
            
            
            figure(h_err3)
            subplot(221)
            plot(e_P_CI_root(:,i_agent)); hold on
            title('trace(P_{cen})/trace(P_{bar})')
            grid on; grid minor;
            subplot(222)
            plot(e_P_bar_root(:,i_agent)); hold on
            title('trace(P_{cen})/trace(P_{CI})')
            grid on; grid minor;
            subplot(2,2,[3 4])
            plot(e_P_bar_cen_2(:,i_agent)); hold on
            plot(e_P_bar_CI_cen_2(:,i_agent)); hold on
            grid on; grid minor;
           
%             plot(e_P_bar_CI_cen_2(:,i_agent)); hold on
            
            
        end
        h_comp =figure;
        for i_agent = 1:9
            figure(h_comp)
            if i_agent==5
                plot(e_bar_BC_dist(:,i_agent),'b','LineWidth',2); hold on
            elseif i_agent==8
                plot(e_bar_BC_dist(:,i_agent),'--r','LineWidth',2); hold on
            end
        end
        xlabel('step')
        ylabel('D_{B}(p_{cen},p_{i})')
         title('Bhatacharia Distance')
         legend('recpetor no. 5','recpetor no. 8')
         grid on
        
        mean(con_perc_ci(range_))
        mean(con_perc_bar(range_))
        mean(con_perc_cen(range_))
        
        
    case '4'
        [flag_coon,G] = generate_graph_for_diag(idx_reg,n_degree_graph,prob_of_link_fail)
        
        
    case '5'
        clear all
        close all
        clc
        flag_noise_model = 'absolute'; % ('absolute' | 'relative')
        global opt_dist
        [A,x0,B,C] = create_sys_atmosphere();
        scale = 1;
        x0 = scale*1000*x0;
        B = scale*1000*B;
        problem_def_diag1(A,x0);
        rank(obsv(A,C))
        %         rank(obsv(A,C(1:7,:)))
        x_prior = x0;
        x_prior_single = x0;
        x_gt = x0;
        q_var =10/(scale)^2;
        r_var = 5/(scale)^2;
        Q = q_var*eye(size(A));
        Q_input = q_var*eye(size(B,2));
        p_prior = Q;
        p_prior_single = 100*Q;
        C_single = C(1:8,:);
        h_predict = figure;
        h_estimate = figure;
        
        for i=1:400
            delta_noise =  randn(size(x_prior))*sqrt(q_var);
            x_predict = A*x_prior +  B*( opt_dist.source.Q');
            x_predict_single = A*x_prior_single +  B*( opt_dist.source.Q');
            P_predict = A*p_prior*A' + Q;
            %             P_predict = A*p_prior*A' + B*Q_input*B';
            P_predict_single = A*p_prior*A' + Q;
            rcond_predcit(i) = rcond(P_predict);
            rcond_predcit_single(i) = rcond(P_predict_single);
            
            disp(['step = ',num2str(i),'---rcond(P_predict) = ',num2str(rcond(P_predict)),'---rcond(P_predict_single)) = ',num2str(rcond(P_predict_single))])
            
            x_gt = A*x_gt +  B*( opt_dist.source.Q')+ delta_noise;
            %             x_gt = A*x_gt +  B*( opt_dist.source.Q' + sqrt(q_var) );
            
            
            %       R = 0.1.*(C*x_gt);
            %             R = r_var*eye(size(C,1));
            %             R_single = r_var*eye(size(C_single,1));
            switch flag_noise_model
                case 'absolute'
                    R = r_var*ones(size(C,1),1);
                    R_single =r_var*ones(size(C_single,1),1);
                case 'relative'
                    R = (0.1.*C*x_gt).^2;
                    R_single = (0.1.*C_single*x_gt).^2;
            end
            %             z = C*x_gt + randn(size(C*x_gt))*sqrt(r_var);
            %             z = C*x_gt + randn(size(C*x_gt)).*sqrt(r_var);
            
            
            z = C*x_gt + randn(size(C*x_gt)).*R;
            h = C*x_predict;
            H = C;
            
            %             z_single = C_single*x_gt + randn(size(C_single*x_gt))*sqrt(r_var);
            z_single = C_single*x_gt + randn(size(C_single*x_gt)).*sqrt(R_single);
            
            %       z = C*x_gt + randn(size(C*x_gt)).*R;
            h_single = C_single*x_predict_single;
            H_single = C_single;
            
            
            %       K = P_predict*H'*inv(H*P_predict*H' + diag(R));
            %             K = P_predict*H'*inv(H*P_predict*H' + (R));
            %             K_single = P_predict_single*H_single'*inv(H_single*P_predict_single*H_single' + (R_single));
            %
            
            K = P_predict*H'*inv(H*P_predict*H' + diag(R));
            K_single = P_predict_single*H_single'*inv(H_single*P_predict_single*H_single' + diag(R_single));
            
            x_est  = x_predict + K*(z-h);
            x_est_single  = x_predict_single + K_single*(z_single-h_single);
            
            p_est = (eye(size(P_predict)) - K*H)*P_predict;
            p_est_single = (eye(size(P_predict_single)) - K_single*H_single)*P_predict_single;
            p_est = 0.5.*(p_est+p_est');
            p_est_single = 0.5.*(p_est_single+p_est_single');
            rcond_est(i) = rcond(p_est);
            rcond_est_single(i) = rcond(p_est_single);
            p_error(i) = norm(p_est - P_predict);
            p_error_single(i) = norm(p_est_single - P_predict_single);
            
            %             err_(i) = norm(x_est -x_gt );
            %             err_single(i) = norm(x_est_single -x_gt );
            
            err_(i) = sqrt(immse(x_est,x_gt)) ;
            err_single(i) = sqrt(immse(x_est_single,x_gt));
            erro_vec(:,i) = x_est -x_gt;
            erro_single_vec (:,i)= x_est_single -x_gt;
            
            trace_p_cen(:,i) = diag(p_est);
            trace_p_single(:,i) = diag(p_est_single);
            
            disp(['step = ',num2str(i),'---rcond(p_est) = ',num2str(rcond(p_est))])
            disp(['step = ',num2str(i),'---norm(p_est - p_prior) = ',num2str(norm(p_est - p_prior))])
            
            disp(['step = ',num2str(i),'---rcond(p_est_single) = ',num2str(rcond(p_est_single))])
            disp(['step = ',num2str(i),'---norm(p_est_single - p_prior_single) = ',num2str(norm(p_est_single - p_prior_single))])
            %       figure(h_predict)
            %       subplot(121); imagesc(P_predict); hold on;
            %       subplot(122); imagesc(p_est); hold on;
            %       pause(0.1)
            %             hold off
            
            x_prior = x_est;
            p_prior = p_est;
            x_prior_single = x_est_single;
            p_prior_single = p_est_single;
            
            perf_index(i) = trace(p_est)/trace(p_est_single);
            bc_dist(i) = BC_distance(x_gt,p_est,x_gt,p_est_single);
            
            
        end
        figure
        subplot(311)
        plot([100:400],rcond_est(100:end)); hold on
        plot([100:400],rcond_est_single(100:end)); hold on
        legend('reciprocal condition number of P_{est}','reciprocal condition number of P_{est_{single}}')
        subplot(312)
        plot(p_error); hold on;
        plot(p_error_single)
        subplot(313)
        plot(err_); hold on;
        plot(err_single); hold on;
        legend('||e_{cen}||','||e_{single}||')
        
        
        figure
        subplot(211)
        plot(perf_index); hold on
        title('performance index centralized vs. single agent')
        
        subplot(212)
        plot(bc_dist); hold on
        title('BC distance bw. centralized and single agent')
        
        h1 = figure;
        h1.Name = 'centralized'
        
        for i_ag=1:9
            subplot(9,1,i_ag)
            plot(erro_vec(i_ag,:),'b'); hold on;
            plot(2.*sqrt(trace_p_cen(i_ag,:)) ,'r'); hold on;
            plot(-2.*sqrt(trace_p_cen(i_ag,:)) ,'r'); hold on;
            
        end
        h2 = figure;
        h2.Name = 'single'
        %         subplot(212)
        for i_ag=1:9
            subplot(9,1,i_ag)
            
            plot(erro_single_vec(i_ag,:),'b'); hold on;
            plot(2.*sqrt(trace_p_single(i_ag,:)) ,'r'); hold on;
            plot(-2.*sqrt(trace_p_single(i_ag,:)) ,'r'); hold on;
        end
        
        
end


