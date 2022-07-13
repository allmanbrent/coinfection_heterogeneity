function void = run_N_script(params)

% This is the script that calls the functions that make the population
% evolve through each time step

gene_mut = zeros(params.N, params.n_genes); % each GENE is initialized with zero mutations
offsp_vec = zeros(params.N);

% N x n_genes x t matrices the store our variables of interest at every time step
num_mut_mat = NaN*zeros(params.N, params.n_genes, params.tstep); % the num of mutations at each time step GENE
fitness_mat = NaN*zeros(params.N, params.n_genes, params.tstep); % the GENE fits at each time step 
virion_fitness_mat = NaN*zeros(params.N, params.tstep); % (N x t) the virion fitnesses at each time step
Npop_mat = NaN*zeros(params.N, params.tstep); % will be the cell number to which each of the virions belong
cell_fits_mat = NaN*zeros(params.C, params.tstep); %the fitness of the cells
offspring_variance = NaN*zeros(1, params.tstep);

for t = 1:params.tstep
    %t
    %log the mutation data from the previous (or initial) generation
    num_mut_mat(:,:,t) = gene_mut; % the mutations on each gene of free virus
    gene_fit = (1-params.sd).^(gene_mut.^(1-params.epistasis));
    fitness_mat(:,:,t) = gene_fit; % the mutations on each gene of free virus
    virion_fitness_mat(:,t) = prod(gene_fit,2);
    
    % distribute virions across cells, and count up how many per cell
    Npop = randi(params.C, params.N, 1); % Npop is a vector of length N that keeps the cell number to which each of the virions belong
    Npop_mat(:,t) = Npop;
    
    % moved mutations and ablation into cells, since seems like they should happen there
    % for ablation, it doesn't matter in terms of output, but seems to me that easier to interpret that way
    %[prog_gene_mut, prog_gene_fit] = mutate(params, gene_mut); %outputs the progeny gene fitnesses and mutation numbers; % this looks good - but mutation occurs to each virion (regardless of whether entered cell or not, and prior to any kind of replication)
    %[abl_mut, abl_fit] = ablateSegments(params); %determine what segments are lost
    
    %do selection with reassortnemt
    [cell_fitness, gene_mut, offsp_vec] = doSelection(t, params, gene_mut, Npop);
    cell_fits_mat(:,t) = cell_fitness;
    offspring_variance(t) = var(offsp_vec);
    
    if max(cell_fitness) == 0
        disp(strcat('No productively infected cells anywhere at time = ', int2str(t)))
        disp('exiting')
        break;
    end
end

% save the simulation results
save(params.outfile, 'params', 'num_mut_mat', 'fitness_mat', 'virion_fitness_mat', 'Npop_mat', 'cell_fits_mat', 'offspring_variance', '-v7.3');
