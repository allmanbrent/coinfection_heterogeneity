function void = GIF_test(void)

%PSIP = 0;
%alpha = 1000;
C=1000;
N=5995;
for r = 1:5
    
    NSEGS=1;%1;
    NGENES=8;
    
    %infile = strcat('alpha_', num2str((alpha)), '_outfile_nsegs', int2str(NSEGS), '_ngenes', int2str(NGENES), '_N', int2str(N), '_C', int2str(C),'_pSIP', (num2str(PSIP*100)), '_r', int2str(r), '.mat')
    infile = strcat('log_NaN_outfile_nsegs1_ngenes8_N',int2str(N),'_C',int2str(C),'_pSIP0_r',int2str(r),'.mat');
    filename = strcat('log_max_n0_nsegs1_ngenes8_N',int2str(N),'_C',int2str(C),'_pSIP0_r',int2str(r),'.gif');
    load(infile);
    if params.U == 0
        filename = strcat('log_max_noMut_var_nsegs1_ngenes8_N',int2str(N),'_C',int2str(C),'_pSIP0_r',int2str(r),'.gif');
    end
    num_mut = reshape((sum(num_mut_mat, 2)), N, params.tstep);
    %hold on; %histogram(num_mut(:,100)); histogram(num_mut(:,200)); histogram(num_mut(:,300)); histogram(num_mut(:,400)); histogram(num_mut(:,500))
    MOI = params.N / params.C;
    %h = figure;
    myplot = figure;
     subplot(2,1,1)
         yyaxis left
         plot(1:params.tstep, mean(num_mut,1));
         title(strcat('MOI=',num2str(MOI),'; N=',num2str(params.N),'; r=', num2str(params.ran)));
         xlabel('time'); ylabel('mean number of mutations');
%     
         yyaxis right
         least_loaded = NaN * zeros(params.tstep, 1);
         for i = 1:params.tstep
             minimum = min(num_mut(:, i));
             least_loaded(i) = length(find(num_mut(:,i) == minimum));
         end
         plot(1:params.tstep, log(least_loaded/params.N));
         ylabel('log fraction least-loaded');
        %plot(1:150, offspring_variance, '--');
        %ylabel('variance in offspring distn');
        %plot(1:150, skewness(num_mut, 1,1), '--');
        %ylabel('skewness of mutation distn');
        %plot(1:150, var(num_mut,1), '--');
        %ylabel('variance of mutation distn');
     subplot(2,1,2)
    %
    
    % PLOT THE HISTOGRAM of MUTATIONS OVER TIME
    %axis tight manual % this ensures that getframe() returns a consistent size
    %'MOI3_nsegs1_pSIP0_r15_fine.gif';
    %filename =strcat('alpha_', num2str(alpha),'_nsegs',num2str(NSEGS),'_ngenes',num2str(NGENES),'_N', int2str(N),'_C',int2str(C),'pSIP',num2str(PSIP*100),'_r',num2str(r),'_mut_hist.gif');
    for t = 1:params.tstep
        % Draw plot for y = x.^n
        
        histogram(num_mut(:,t), [0:1:200])%,'NumBins', 15);
        if (MOI >= 1) && (MOI < 5)
            histogram(num_mut(:,t), [0:1:200])
            axis([0 180 0 N/10]) %
        elseif MOI < 1
            histogram(num_mut(:,t), [0:1:50])
            axis([0 60 0 450])
        else
            histogram(num_mut(:,t), [0:1:300])
            axis([0 300 0 N/10])
        end
        xlabel('number of mutations in genome'); ylabel('counts of individuals');
        %text(35, 100, sprintf(strcat('t = ', int2str(t))));
        %title(strcat('t= ', int2str(t), ' MOI=',num2str(MOI),'; ', int2str(NSEGS),' segment(s);', int2str(NGENES), ' genes; r=', num2str(r), '; alpha=', num2str(alpha)));
        title(strcat('t= ', int2str(t), '; MOI=',num2str(MOI),'; r=', num2str(params.ran)));
        %plot(x,y)
        xline(min(num_mut(:,t)));
        xline(mean(num_mut(:,t)), 'r');
        drawnow
        % Capture the plot as an image
        %frame = getframe(h);
        frame = getframe(myplot(:));
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        % Write to the GIF File
        if t == 1
            imwrite(imind,cm,filename,'gif','DelayTime', 0.5, 'Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif', 'DelayTime', 0.1,'WriteMode','append');
        end
    end
end

N = 21544;
C = 1000;
myplot2 = figure;
hold on;
for ran = 6:8
    curr_file = strcat('log_NaN_outfile_nsegs1_ngenes8_N', num2str(N),'_C', num2str(C),'_pSIP0_r', num2str(ran), '.mat');
    load(curr_file);
    num_mut = reshape((sum(num_mut_mat, 2)), N, params.tstep);
    plot(1:150, mean(num_mut,1));
    title(strcat('t= ', int2str(t), ' MOI=',num2str(N/C),'; N=',num2str(params.N)));
    xlabel('time'); ylabel('mean number of mutations')
end
hold off

end