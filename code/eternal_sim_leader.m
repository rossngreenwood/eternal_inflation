
function eternal_sim_leader(cores,flag_range,range,input_file,output_dir)
    
    warning('off');
    
    if nargin < 4, input_file = 'infile.txt';   end
    if nargin < 5, output_dir = './';           end
    
    parpool(cores) % Open a parallel pool
    
    if flag_range
        test_ids = range(1):range(2);
    else
        test_ids = range;
    end
    
    for test_id = test_ids
        
        %% Gather parameters from input file
        
        fid = fopen(sprintf(input_file,test_id),'r');
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
        [mv_offset_max,~,~,is1]     = sscanf(meta_line(is:end),'%G,',1); is = is+is1-1;
        [aux_flag,~,~,is1]          = sscanf(meta_line(is:end),'%d,',1); is = is+is1-1;
        [Nafter,~,~,is1]            = sscanf(meta_line(is:end),'%G,',1); is = is+is1-1;
        [seed,~,~,is1]              = sscanf(meta_line(is:end),'%d,',1); is = is+is1-1;
        [n_recycle,~,~,is1]         = sscanf(meta_line(is:end),'%d,',1); is = is+is1-1;
        
        if seed == -1
            rs = RandStream.create('mrg32k3a','NumStreams',cores,'Seed','shuffle','CellOutput',true);
        else
            rs = RandStream.create('mrg32k3a','NumStreams',cores,'Seed',seed,'CellOutput',true);
        end
        
        %% Run code in parallel
        
        spmd % Execute this code in parallel
            
            worker_outfile = strcat(sprintf(output_dir,test_id),...
                sprintf('.worker_%d.txt',labindex-1));
            
            warning('off');
            
            eis = EternalInflationSimulator(...
                'outfile',          worker_outfile,...
                'n_iter',           floor(n_iter/cores),...
                'mv',               mv,...
                'mh',               mh,...
                'kmax',             kmax,...
                'gamma',            gamma,...
                'measure',          measure,...
                'n_tunnel_max',     n_tunnel_max,...
                'lambdascreenmode', logical(lambdascreen),...
                'mv_offset_max',    mv_offset_max,...
                'Nafter',           Nafter,...
                'n_recycle',        n_recycle,...
                'randstream',       rs{labindex},...
                'cores',            cores,...
                'flag_screen_basin',aux_flag);
            
            eis.main();
            
        end
        
    end % for

    delete(gcp('nocreate'));
    
end
