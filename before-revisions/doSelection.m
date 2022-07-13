function [cell_fitness, next_gene_mut, offsp_vec] = doSelection(t, params, gene_mut, Npop)

% moved the ablation step here first, prior to cell fitness calculation
cell_fitness = NaN*zeros(params.C, 1);
infection_sz = NaN*zeros(params.C, 1);
for i = 1:params.C
    cell_virus(i).locs = find(Npop == i);
    cell_virus(i).n_virions = length(cell_virus(i).locs);
    
    %we store this number so we can scale the output later
    infection_sz(i) = length(cell_virus(i).locs);
    
    % copy gene mut
    cell_virus(i).gene_mut = gene_mut(cell_virus(i).locs, :);
    % Ablate the segments in this cell
    [cellular_segment_ablation, cellular_gene_mut] = AblateSegments(params, cell_virus(i).n_virions, cell_virus(i).gene_mut);
    cell_virus(i).segment_ablation = cellular_segment_ablation; % the ablated segments
    cell_virus(i).gene_mut = cellular_gene_mut; % the genes with ablation
    
    cell_fitness(i,1) = GetCellFitness(params, cell_virus(i));
end

if max(cell_fitness) == 0
    disp(strcat('No productively infected cells anywhere at time = ', int2str(t)))
    disp('exiting')
    gene_mut = []; gene_fit  = [];
    return;
end

if params.heterogeneous == true
    %successful_cell_indices = dirmultinom(cell_fitness, params, t);
    successful_cell_indices = gammaHet(cell_fitness, params);
elseif params.scale_out_lin == true
    successful_cell_indices = outLinear(cell_fitness, infection_sz, params);
%elseif params.scale_out_exp == true
 %   successful_cell_indices = outExponential(cell_fitness, infection_sz, params);
else
    rel_fitness = cell_fitness/sum(cell_fitness); %relativize cellular fitness
    % perform a multinomial sampling process:
    successful_cell_indices = randsample(1:params.C, params.N, true, rel_fitness); % N draws from the vals 1 to C, based on relative fitnesses
end 

out_per_cell = tabulate(successful_cell_indices);
out_per_cell = out_per_cell(:,2);
cell_output_counts = cell(params.C, params.n_segments);

% then mutate at the step of generating output virus
for j = 1:params.N %create N progeny from the multinomial draws
    cell_index = successful_cell_indices(j);
    [next_gene_mut(j,:), indivs_chosen] = GetOutputVirus(cell_virus(cell_index), params);
    if params.n_segments == 1
        account = tabulate(indivs_chosen); %the numbers of indivs chosen of each within cell index
        account = [account(:,2)' zeros(1, infection_sz(cell_index) - max(indivs_chosen))];
        if size(cell_output_counts{cell_index}) == [0 0]
            cell_output_counts{cell_index} = account;
        else
            cell_output_counts{cell_index} = cell_output_counts{cell_index} + account;
        end
    end
end

%the offspring distribution making sure to fill the rest of the cells such
%that all of the virions that did not reproduce and are singly infected
%also represented
if params.n_segments == 1
    offsp_vec = [[cell_output_counts{:}] zeros(1,params.N - length([cell_output_counts{:}]))];
else
    offsp_vec = NaN;
end
