function void = params_segmented_figs()

%%Figure 4 SEGMENTATION BASE CASE
    F4params.heterogeneous = false;
    F4params.scale_out_lin = false;
    F4params.k = 0;
    F4params.ran = 1:20;
    F4params.U = 1;
    F4params.epistasis = 0;
    F4params.sd = 0.2;
    F4params.tstep = 150;
    F4params.n_segments = [2 4 8];
    F4params.pSIP = 0;
    F4params.n_genes = 8;

    %%Figure 4A
    F4params.N = 1000;
    MOI = 0.1;
    F4params.C = round(F4params.N ./ MOI);
    main_coinfection_script(F4params);
    summarize_data_NeV(F4params);
    
    %%Figure 4B
    F4params.N = round(logspace(log10(100),log10(10000), 10));
    MOI = 0.1;
    for n = F4params.N
        F4params.C = round(n ./ MOI);
        for ns = F4params.n_segments
            params = F4params;
            params.n_segments = ns;
            main_coinfection_script(params);
            summarize_data_NeV(params);
        end
    end
    
    %%Figure 4C
     F4params.N = 1000;
     MOI = logspace(log10(.01),log10(1000), 10);
     F4params.C = round(F4params.N ./ MOI);
%     main_coinfection_script(F4params);
     for c = F4params.C
         for n = F4params.n_segments
             params = F4params;
             params.C = c;
             params.n_segments = n;
             main_coinfection_script(params);
             summarize_data_NeV(params);
         end
     end     
     
    %%Figure 4D
    F4params.C = 1000;
    MOI = logspace(log10(.01),log10(1000), 10);
    F4params.N = round(MOI .* F4params.C);
    for n = F4params.N
        for ns = F4params.n_segments
            params = F4params;
            params.N = n;
            params.n_segments = ns;
            main_coinfection_script(params);
            summarize_data_NeV(params);
        end
    end

end