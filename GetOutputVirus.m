function [next_virus_by_gene, indivs_chosen] = GetOutputVirus(cell_info, params)

indivs_chosen = NaN*zeros(params.n_segments);
% do reassortment in here
for i = 1:params.n_segments
    locsNotNaN = find(~isnan(cell_info.segment_ablation(:,i)));
    %pick a number from 1 to # indivs w/o ablated segment i
    loc_index = randi(length(locsNotNaN)); 
    locChosen = locsNotNaN(loc_index); %choose a segment id that wasn't ablated
    genes = find(params.gene_on_segment_vector == i);
    %pick the genes on that chosen segment
    next_virus_by_gene(1, genes) = cell_info.gene_mut(locChosen, genes);
    indivs_chosen(i) = locChosen;
end


% then mutate

perGeneU = params.U/params.n_genes; % the per gene mutation rate
% how many mutations to add to each gene
newMut_matrix = poissrnd(perGeneU, size(next_virus_by_gene)); 
% add the new mutations to the parental mutation number
next_virus_by_gene = next_virus_by_gene + newMut_matrix; 
