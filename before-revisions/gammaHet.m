function [successful_cell_indices] = gammaHet(cell_fits, params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
rel_cell_fits = cell_fits/sum(cell_fits);
%must calculate theta as the shape parameter for the gamma distn
theta = rel_cell_fits./params.k; 

sum_gamma_rands = 0;
numIterations = 0;
maxIterations = 1000; % More than you ever expect to have.
while (sum_gamma_rands == 0) && (numIterations < maxIterations)
   gammarands = gamrnd(params.k, theta); %they use the k,theta 
   % Increment loop counter.
   numIterations = numIterations + 1;
   % Recompute condition for continuing.
   sum_gamma_rands = sum(gammarands); %if all of the numbers are zero, choose a new set
end

%relativize the gamma numbers, which will be the probabilities in our
%multinomial random draw
multinom_probs = gammarands / sum(gammarands);

% N draws from the vals 1 to C, based on relative fitnesses
successful_cell_indices = randsample(1:params.C, params.N, true, multinom_probs); 

end

