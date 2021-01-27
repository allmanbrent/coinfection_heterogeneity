function [cell_fitness] = GetCellFitness(params, cell_info)

% moved your cell fitness calculations into this function

if cell_info.n_virions == 0 % no virus got in
    cell_fitness = 0;
    return;
end

%identify the ablated genes
n_genes_ablated_in_cell = sum(isnan(cell_info.gene_mut), 1);

if any(n_genes_ablated_in_cell == cell_info.n_virions)
    % does not contain essential set of genes
    cell_fitness = 0;
    return;
end

gene_fitness = (1-params.sd).^(cell_info.gene_mut.^(1-params.epistasis));

% this dimension takes the means for each gene
mean_gene_fits = mean(gene_fitness, 1, 'omitnan'); 

cell_fitness = prod(mean_gene_fits); %the cellular fitness is the geometric mean of the avg segment fitnesses

