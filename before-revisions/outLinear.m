function [successful_cell_indices] = outLinear(cell_fits, infection_sz, params)

%scale the cell fitnesses linearly with the number
    % of virions infecting that cell
scaled_cell_fits = cell_fits .* infection_sz;
rel_cell_fits = scaled_cell_fits/sum(scaled_cell_fits);

%the rel_cell_fits act as multinomial probabilities scaled by the # of
%virions in the cell
successful_cell_indices = randsample(1:params.C, params.N, true, rel_cell_fits); 
end

