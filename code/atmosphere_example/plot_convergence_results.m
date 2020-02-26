function plot_convergence_results(result,txt_)
figure
subplot(121)
global nCovSamples;

for i=1:nCovSamples
    for j=1:70
        if strcmp(txt_,'ICI_trace') || strcmp(txt_,'ICI_det')
            det_(i,j) = det(inv(((result.consenus{i}.P_prior{j}))));
        else
            det_(i,j) = det(inv(((result.consenus{i}.Y_prior{j}))));
        end
        
    end
end

for i =1:nCovSamples
    plot(det_(i,:)); hold  on
end
switch txt_
    case 'ICI_trace'
        title(['Convergence of ','ICI_{trace}'])

    case 'ICI_det'
        title(['Convergence of ','ICI_{det}'])

    case 'MH'
        title(['Convergence of ','MH'])

end        


        
% title(['Convergence of',txt_])
xlabel('steps')
ylabel('det(P)')
grid on
grid minor




subplot(122)
for i=1:nCovSamples
    for j=1:70
        if strcmp(txt_,'ICI_trace') || strcmp(txt_,'ICI_det')
            det_(i,j) = trace(inv(((result.consenus{i}.P_prior{j}))));
        else
            det_(i,j) = trace(inv(((result.consenus{i}.Y_prior{j}))));
        end
    end
end


for i =1:nCovSamples
    plot(det_(i,:)); hold  on
end

switch txt_
    case 'ICI_trace'
        title(['Convergence of ','ICI_{trace}'])

    case 'ICI_det'
        title(['Convergence of ','ICI_{det}'])

    case 'MH'
        title(['Convergence of ','MH'])
end        


% title(['Convergence of',txt_])
xlabel('steps')
ylabel('trace(P)')
grid on
grid minor


