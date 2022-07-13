function void = params_min_model_figs()

%%Figure 1 BASE MODEL
    F1params.heterogeneous = false;
    F1params.scale_out_lin = false;
    F1params.k = 0;
    F1params.ran = 1:10;
    F1params.U = 1;
    F1params.epistasis = 0;
    F1params.sd = 0.2;
    F1params.tstep = 150;
    F1params.n_segments = 1;
    F1params.pSIP = 0;
    F1params.n_genes = 8;

%Figure 1A
%     F1params.N = 1000;
%     F1params.C = 10000;
%     main_coinfection_script(F1params);
%     summarize_data_NeV(F1params);
%     
%Figure 1B
     F1params.N = [1 10 round(logspace(log10(100),log10(10000), 10))];
     MOI = 0.1;
     for n = F1params.N
         F1params.C = round(n ./ MOI);
         params = F1params;
         params.N = n;
         main_coinfection_script(params);
         summarize_data_NeV(params);
     end
    %F1params.ran = 25:34;
%Figure 1C
    F1params.N = 1000;
    MOI = logspace(log10(.01),log10(1000), 10);
    F1params.C = round(F1params.N ./ MOI);
    for n = F1params.C
        params = F1params;
        params.C = n;
        main_coinfection_script(params);
        summarize_data_NeV(params);
    end
    
%Figure 1D
    F1params.C = 1000;
    MOI = logspace(log10(.01),log10(1000), 10);
    F1params.N = round(MOI .* F1params.C);
    for n = F1params.N(10)
        params = F1params;
        params.N = n;
        if n > 77426
            params.ran = 1;
        end
        main_coinfection_script(params);
        summarize_data_NeV(params);
    end
end