function [ablated_segment, ablated_gene_mut] = AblateSegments(params, n_virions, gene_mut)

% ablate only the virions that made it into that cell

abl_mut = zeros(n_virions, params.n_segments);
% create a matrix with random, uniformly distributed values between 0 and 1
rand_vals = rand([n_virions params.n_segments]);
% find the cells that fall below the value of params.pSIP
locsZero = find(rand_vals < params.pSIP);

abl_mut(locsZero) = NaN; %ablate the gene segments
ablated_segment = abl_mut;

ablated_gene_mut = gene_mut;
%determine which genes were lost on the ablated segments
for i = 1:n_virions % loop through the virions
    for j = 1:params.n_segments %loop through the segments
       if isnan(abl_mut(i,j)) % if the segments was lost, lose the genes
           genes_on_segment = find(params.gene_on_segment_vector == j);
           ablated_gene_mut(i, genes_on_segment) = NaN;
       end
    end
end
