function void = params_dep_output_figs()

%%INPUT-DEPENDENT OUTPUT
    F3params.heterogeneous = false;
    F3params.scale_out_lin = true;
    F3params.k = 0;
    F3params.ran = 1:20;
    F3params.U = 1;
    F3params.epistasis = 0;
    F3params.sd = 0.02;
    F3params.tstep = 150;
    F3params.n_segments = 1;
    F3params.pSIP = 0;
    F3params.n_genes = 8;

    %%Figure 3A
     F3params.N = logspace(log10(100),log10(10000), 10);
     MOI = 0.1;
     for n = F3params.N
         params = F3params;
         params.C = round(n ./ MOI);
         params.N = n;
         main_coinfection_script(params);
         summarize_data_NeV(params);
     end

    %%Figure 3B
    F3params.N = 1000;
    MOI = logspace(log10(.01),log10(1000), 10);
    F3params.C = round(F3params.N ./ MOI);
    for n = F3params.C
        params = F3params;
        params.C = n;
        main_coinfection_script(params)
        summarize_data_NeV(params);
    end
    
    %%Figure 3C
    F3params.C = 1000;
    MOI = logspace(log10(.01),log10(1000), 10);
    F3params.N = round(MOI .* F3params.C);
    for n = F3params.N
        params = F3params;
        params.N = n;
        main_coinfection_script(params)
        summarize_data_NeV(params);
    end
end