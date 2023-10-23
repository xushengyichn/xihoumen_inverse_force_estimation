load results_experiment.mat

clear t u udot u2dot p_filt_m h_hat udot_h power total_work_done cumulated_work_done

for t1 = 1:length(experiment_names)
    
    exp_name = experiment_names{t1};
    t(:,t1) = results_experiment.(exp_name).t;
    u(:,t1) = results_experiment.(exp_name).u;
    udot(:,t1) = results_experiment.(exp_name).udot;
    u2dot(:,t1) = results_experiment.(exp_name).u2dot;
    p_filt_m(:,t1) = results_experiment.(exp_name).p_filt_m;
    h_hat = results_experiment.(exp_name).h_hat;
    udot_h(:,t1)=h_hat(3,:)/2.008345587515120e-04;
    t_seconds = seconds(t(:,t1) - t(1,t1));
    % power(:,t1)=p_filt_m(:,t1).*h_hat(3, :)';
    power(:,t1)=p_filt_m(:,t1).*udot_h(:,t1);
    total_work_done(t1) = trapz(t_seconds, power(:,t1));
    cumulated_work_done(:,t1) = cumtrapz(t_seconds, power(:,t1));
end

figure
plot(t(:,1),cumulated_work_done(:,1))
hold on
plot(t(:,2),cumulated_work_done(:,2))
plot(t(:,3),cumulated_work_done(:,3))
plot(t(:,4),cumulated_work_done(:,4))

figure
plot(t(:,1),power(:,1))
hold on
plot(t(:,2),power(:,2))
plot(t(:,3),power(:,3))
plot(t(:,4),power(:,4))