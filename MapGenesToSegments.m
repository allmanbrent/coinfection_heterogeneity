function gene_on_segment_vector = MapGenesToSegments(params)

g_per_s = params.n_genes/params.n_segments; % the number of genes per gene segment

for i = 1:params.n_segments
    start_gene = 1+(i-1)*g_per_s;
    end_gene = i*g_per_s; 
    gene_on_segment_vector(1, start_gene:end_gene) = i;
end            