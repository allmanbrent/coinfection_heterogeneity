function void = summarize_data_NeV(input_params)

    params = struct();
    params.heterogeneous = input_params.heterogeneous;
    params.scale_out_lin = input_params.scale_out_lin;
    params.tstep = input_params.tstep;

for n = 1:length(input_params.N)

params.N = input_params.N(n);
    
for c = 1:length(input_params.C)
    if length(input_params.C) == 1
        params.C = input_params.C(1);
    elseif length(input_params.C) == length(input_params.N) && length(input_params.N) > 1
        %we don't want to loop over all C values for all N values
        params.C = input_params.C(n); 
    else
        params.C = input_params.C(c);
    end
for NSEGS = input_params.n_segments
    params.n_segments = NSEGS;
for NGENES = input_params.n_genes
    params.n_genes = NGENES;
    % do this once here, instead of many
    % times in reassortment function
    params.gene_on_segment_vector = MapGenesToSegments(params);
for E = input_params.epistasis
    params.epistasis = E;
for SD = input_params.sd
    params.sd = SD;
for PSIP = input_params.pSIP
    params.pSIP = PSIP;
    params.proportionVirusSIP = 1 - (1-params.pSIP)^params.n_segments;
for mu = input_params.U
    params.U = mu;
for K = input_params.k
    params.k = K;

    mean_per_t_per_rep = zeros(20, max(input_params.ran));%zeros(input_params.tstep, max(input_params.ran));
    offspring_var_means = NaN(1, max(input_params.ran));
    mean_per_t = zeros(150);%zeros(params.tstep);
    for RAN = input_params.ran
        params.ran = RAN;
        if params.epistasis == 0 && params.heterogeneous == true && params.scale_out_lin == false
            params.outfile = strcat('gamma_', num2str(round(params.k*100)), '_outfile_nsegs', int2str(params.n_segments), '_ngenes', int2str(params.n_genes), '_N', int2str(params.N), '_C', int2str(params.C),'_pSIP', (num2str(params.pSIP*100)), '_r', int2str(params.ran))
            load(params.outfile);
            if params.N == 1
                mean_per_t_per_rep(:,RAN) = reshape((sum(num_mut_mat, 2)), params.tstep,1);
            else
                mean_per_t_per_rep(:,RAN) = reshape(mean(sum(num_mut_mat(:,:,1:20), 2)), 20,1);
                if exist('offspring_variance') == 1
                    offspring_var_means(RAN) = mean(offspring_variance);
                end
            end
            outfile_name = strcat('summary_gamma_', num2str(round(params.k*100)), '_nsegs', int2str(params.n_segments), '_ngenes', int2str(params.n_genes), '_N', int2str(params.N), '_C', int2str(params.C),'_pSIP', (num2str(params.pSIP*100)));
        elseif params.epistasis == 0 && params.heterogeneous == false && params.scale_out_lin == false
            params.outfile = strcat('log_NaN_outfile_nsegs', int2str(params.n_segments), '_ngenes', int2str(params.n_genes), '_N', int2str(params.N), '_C', int2str(params.C),'_pSIP', (num2str(params.pSIP*100)), '_r', int2str(params.ran))
            load(params.outfile);
            if params.n_segments > 1
                mean_per_t_per_rep(:,RAN) = reshape(mean(sum(num_mut_mat, 2)), params.tstep,1);
            else
                mean_per_t_per_rep(:,RAN) = reshape(mean(sum(num_mut_mat, 2)), params.tstep,1);
            end
            outfile_name = strcat('summary_log_NaN_nsegs', int2str(params.n_segments), '_ngenes', int2str(params.n_genes), '_N', int2str(params.N), '_C', int2str(params.C),'_pSIP', (num2str(params.pSIP*100)));
        elseif params.scale_out_lin == true
            params.outfile = strcat('lin_outfile_nsegs', int2str(params.n_segments), '_ngenes', int2str(params.n_genes), '_N', int2str(params.N), '_C', int2str(params.C),'_pSIP', (num2str(params.pSIP*100)), '_r', int2str(params.ran))
            load(params.outfile);
            mean_per_t_per_rep(:,RAN) = reshape(mean(sum(num_mut_mat, 2)), params.tstep,1);
            outfile_name = strcat('summary_lin_nsegs', int2str(params.n_segments), '_ngenes', int2str(params.n_genes), '_N', int2str(params.N), '_C', int2str(params.C),'_pSIP', (num2str(params.pSIP*100)));
        else
            params.outfile = strcat('epsilon_', int2str(params.epistasis*100), '_nsegs', int2str(params.n_segments), '_ngenes', int2str(params.n_genes), '_N', int2str(params.N), '_C', int2str(params.C),'_pSIP', (num2str(params.pSIP*100)), '_r', int2str(params.ran))
            load(params.outfile);
            mean_per_t_per_rep(:,RAN) = reshape(mean(sum(num_mut_mat, 2)), params.tstep,1);
            
        end
    end
mean_per_t = mean(mean_per_t_per_rep, 2);
std_error = std(mean_per_t_per_rep,0,2) / sqrt(max(input_params.ran));

if isnan(offspring_var_means(1)) == 0
    V_eff_mean = mean(params.N ./ offspring_var_means);
    V_eff_error = std(offspring_var_means) ./ sqrt(max(input_params.ran));
    save(outfile_name, 'params', 'mean_per_t', 'std_error', 'V_eff_mean', 'V_eff_error');
else
    save(outfile_name, 'params', 'mean_per_t', 'std_error');
end

end
end
end
end
end
end
end
end

end
end
