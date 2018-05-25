
function eternal_sim_leader(cores,input_file,output_file,output_dir)
    
    warning('off');
    
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
    
    if seed == -1
        rs = RandStream.create('mrg32k3a','NumStreams',cores,'Seed','shuffle','CellOutput',true);
    else
        rs = RandStream.create('mrg32k3a','NumStreams',cores,'Seed',seed,'CellOutput',true);
    end
    
    %% Run code in parallel
    
    parpool(cores) % Open a parallel pool
    
    spmd % Execute this code in parallel
        
        worker_outfile = strcat(output_dir,sprintf('.worker_%d.txt',labindex-1));
        
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
            'rho_Lambda_thres', rho_Lambda_thres,...
            'fixQ',             logical(fixQ),...
            'Nafter',           Nafter,...
            'n_recycle',        n_recycle,...
            'randstream',       rs{labindex});
        
        eis.main();
        
    end
    
end
