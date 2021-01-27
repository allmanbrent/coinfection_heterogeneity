function void = main_coinfection_script(input_params)

%clear all; close all; clc;

% params.heterogeneous = input_params.heterogeneous;
% params.scale_out_lin = input_params.scale_out_lin;
% params.tstep = input_params.tstep;

parfor RAN = input_params.ran
    temp_ran = RAN;
    params = struct();
    params.heterogeneous = input_params.heterogeneous;
    params.scale_out_lin = input_params.scale_out_lin;
    params.tstep = input_params.tstep;
    params.ran = temp_ran;

for n = input_params.N
%     if params.ran == 9 || params.ran == 17
%         if n >= 5595
%             params.N = n;
%         else
%             continue
%         end
%     else
%         params.N = n;
%     end
params.N = n;
for c = input_params.C
    params.C = c;
for NSEGS = input_params.n_segments;
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

if params.epistasis == 0 && params.heterogeneous == true && params.scale_out_lin == false
    %for logmean = LOGMEAN
    %for logvar = LOGVAR%ALPHA = alpha
    %params.lognormean = logmean;
    %params.lognormvar = logvar;%params.alpha = ALPHA;
    %params.outfile = strcat('log_', num2str((params.lognormean)),'_', num2str((params.lognormvar)), '_outfile_nsegs', int2str(params.n_segments), '_ngenes', int2str(params.n_genes), '_N', int2str(params.N), '_C', int2str(params.C),'_pSIP', (num2str(params.pSIP*100)), '_r', int2str(params.ran))
    params.outfile = strcat('gamma_', num2str(round(params.k*100)), '_outfile_nsegs', int2str(params.n_segments), '_ngenes', int2str(params.n_genes), '_N', int2str(params.N), '_C', int2str(params.C),'_pSIP', (num2str(params.pSIP*100)), '_r', int2str(params.ran))
    %if isfile(strcat(params.outfile, '.mat')) == 0
        rng(params.ran);
        run_N_script(params); %run the simulation WITH reassortment
    %end
elseif params.epistasis == 0 && params.heterogeneous == false && params.scale_out_lin == false
    params.alpha = NaN;
    params.outfile = strcat('log_NaN_outfile_nsegs', int2str(params.n_segments), '_ngenes', int2str(params.n_genes), '_N', int2str(params.N), '_C', int2str(params.C),'_pSIP', (num2str(params.pSIP*100)), '_r', int2str(params.ran))
    if isfile(strcat(params.outfile, '.mat')) == 0
        rng(params.ran);
        run_N_script(params); %run the simulation WITH reassortment
    end
elseif params.scale_out_lin == true
    params.outfile = strcat('lin_outfile_nsegs', int2str(params.n_segments), '_ngenes', int2str(params.n_genes), '_N', int2str(params.N), '_C', int2str(params.C),'_pSIP', (num2str(params.pSIP*100)), '_r', int2str(params.ran))
    %if isfile(strcat(params.outfile, '.mat')) == 0
        rng(params.ran);
        run_N_script(params); %run the simulation WITH reassortment
    %end
elseif params.scale_out_exp == true
    params.outfile = strcat('exp_outfile_nsegs', int2str(params.n_segments), '_ngenes', int2str(params.n_genes), '_N', int2str(params.N), '_C', int2str(params.C),'_pSIP', (num2str(params.pSIP*100)), '_r', int2str(params.ran))
    if isfile(strcat(params.outfile, '.mat')) == 0
        rng(params.ran);
        run_N_script(params); %run the simulation WITH reassortment
    end
else
    params.outfile = strcat('epsilon_', int2str(params.epistasis*100), '_nsegs', int2str(params.n_segments), '_ngenes', int2str(params.n_genes), '_N', int2str(params.N), '_C', int2str(params.C),'_pSIP', (num2str(params.pSIP*100)), '_r', int2str(params.ran))
    if isfile(strcat(params.outfile, '.mat')) == 0
        rng(params.ran);
        run_N_script(params); %run the simulation WITH reassortment
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
end
end

