function [probs] = verify_output(svms, alpha, data)
    
    plot_this = true;

    num_sets = length(svms);
    probs = zeros(num_sets,1);
    
    Nsamples = size(data,1);

    for i = 1:num_sets
        
        [inside,~] = svms{i}.predict(data);
        probs(i) = sum(inside == 1)/ Nsamples;
    end

    verified = all(probs >= 1 - alpha);

    fprintf('Validity: %d\n', verified);
    %fprintf('alpha values: %f\n', alpha');
    %fprintf('probabilities: %f\n', probs');
    
    if plot_this
        figure()
        stairs(1 - alpha, probs)
        hold on
        plot([0 1], [0 1], 'r')
        xlabel("1 - alphas")
        ylabel("probs")
    end

end