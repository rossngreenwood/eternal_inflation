
function eternal_sim_leader(cores,input_file,output_file,output_dir)
    
    if nargin < 2, input_file = 'infile.txt';   end
    if nargin < 3, output_file = 'outfile.txt'; end
    if nargin < 4, output_dir = '';             end
    
    %% Gather parameters from input file
    
    fid = fopen(input_file,'r');
    meta_line = fgets(fid);
    is = 1; fclose(fid);
    
    [n_iter,~,~,is1] = sscanf(meta_line(is:end),'%E,',1); is = is+is1-1;
    [mv,~,~,is1]     = sscanf(meta_line(is:end),'%G,',1); is = is+is1-1;
    [mh,~,~,is1]     = sscanf(meta_line(is:end),'%G,',1); is = is+is1-1;
    [kmax,~,~,is1]   = sscanf(meta_line(is:end),'%d,',1); is = is+is1-1;
    [gamma,~,~,is1]  = sscanf(meta_line(is:end),'%f,',1); is = is+is1-1;

    [measure,~,~,is1]           = sscanf(meta_line(is:end),'%c,',1); is = is+is1-1;
    [n_tunnel_max,~,~,is1]      = sscanf(meta_line(is:end),'%d,',1); is = is+is1-1;
    [lambdascreen,~,~,is1]      = sscanf(meta_line(is:end),'%d,',1); is = is+is1-1;
    [rho_Lambda_thres,~,~,is1]  = sscanf(meta_line(is:end),'%G,',1); is = is+is1-1;
    [fixQ,~,~,is1]              = sscanf(meta_line(is:end),'%d,',1); is = is+is1-1;
    [Nafter,~,~,is1]            = sscanf(meta_line(is:end),'%G,',1); is = is+is1-1;
    [seed,~,~,is1]              = sscanf(meta_line(is:end),'%d,',1); is = is+is1-1;
    [n_recycle,~,~,is1]         = sscanf(meta_line(is:end),'%d,',1); is = is+is1-1;
    
    disp(n_iter,mv,mh,kmax,gamma,measure,n_tunnel_max,...
        lambdascreen,rho_Lambda_thres,fixQ,Nafter,seed,n_recycle);
    
    eis_wrapper_partial = @(wid) eis_wrapper(wid,...
        worker_iter,...
        mh,...
        mv,...
        kmax,...
        gamma,...
        measure,...
        n_tunnel_max,...
        lambdascreen,...
        rho_Lambda_thres,...
        fixQ,...
        Nafter,...
        seed,...
        n_recycle,...
        output_dir);
        
    parpool(cores)
    spmd
        eis_wrapper_partial(labindex);
    end
    
    disp('Combining output...')
    
%     % Now that all the processes have completed, we need to process the worker specific vcfs
%     worker_files = arrayfun(@(i) sprintf('.worker_%d.txt',i), 0:cores-1,'Un',0);
%     fid = fopen([output_dir output_file],'wt');
%     for filename = worker_files
%         w_fid = fopen([output_dir filename],'r');
%         header_line = w_file.readline();
%         if filename == 0
%             % Only write the header line once
%             print(header_line, file=outfile, end='')
%         end
%         for line in w_file:
%             fprintf(fid,'%s',line)
%         end
%     end
%     for filename in worker_files:
%         % Delete the temp files now that we're done with them.
%         os.remove(params.output_dir + worker_files[filename]);
%     end
    
end

function eis_wrapper(worker_id, worker_iter, mh, mv, kmax, gamma, measure, n_tunnel_max, lambdascreen, rho_Lambda_thres, fixQ, Nafter, seed, n_recycle, output_dir)
    
    worker_id = num2str(worker_id);
    worker_outfile = [output_dir sprintf('.worker_%s.txt',worker_id)];
    
    warning('off');
    
    eis = EternalInflationSimulator(...
        'outfile',worker_outfile,...
        'n_iter',worker_iter,...
        'mv',mv,...
        'mh',mh,...
        'kmax',kmax,...
        'gamma',gamma,...
        'measure',measure,...
        'n_tunnel_max',n_tunnel_max,...
        'lambdascreenmode',logical(lambdascreen),...
        'rho_Lambda_thres',rho_Lambda_thres,...
        'fixQ',logical(fixQ),...
        'Nafter',Nafter,...
        'seed',randi(10000),...
        'n_recycle',n_recycle);
    eis.main();
    
end

