function void = params_stoch_het_figs()

    F2params.heterogeneous = true;
    F2params.scale_out_lin = false;
    F2params.k = logspace(log10(0.01),log10(10000), 5);
    F2params.ran = 1:20;
    F2params.U = 1;
    F2params.epistasis = 0;
    F2params.sd = 0.2;
    F2params.tstep = 20;
    F2params.n_segments = 1;
    F2params.pSIP = 0;
    F2params.n_genes = 8;

    %%Figure 2A
    F2params.N = round(logspace(log10(100),log10(10000), 10));
    MOI = 0.1;
    for n = F2params.N
        F2params.C = round(n ./ MOI);
        params = F2params;
        params.N = n;
        main_coinfection_script(params);
        summarize_data_NeV(params);
    end
    
    %%Figure 2C
    F2params.N = 1000;
    MOI = logspace(log10(.01),log10(1000), 10);
    F2params.C = round(F2params.N ./ MOI);
    for c = F2params.C
        params = F2params;
        params.C = c;
        main_coinfection_script(params);
        summarize_data_NeV(params);
    end
    
    %%Figure 2D
    F2params.C = 1000;
    MOI = logspace(log10(.01),log10(1000), 10);
    F2params.N = round(MOI .* F2params.C);
    for n = F2params.N
        params = F2params;
        params.N = n;
        main_coinfection_script(params);
        summarize_data_NeV(params);
    end
end