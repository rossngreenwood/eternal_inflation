function eis_wrapper(...
        outfile,            ...
        n_iter,             ...   % Number of iterations
        mv,                 ...   % Initial mass scale of the potential (multiple of Mpl)
        mh,                 ...   % Mass scale of the inflaton field (multiple of Mpl)
        kmax,               ...   % Largest wavenumber for GRF
        gamma,              ...   % Frequency dependence of GRF
        measure,            ...   % Measure on initial conditions
        n_tunnel_max,        ...  % Max number of tunneling events to simulate
        lambdascreen,       ...   % Throw out cases where rho_Lambda < 0?
        rho_Lambda_thres,   ...   % Condition on rho_Lambda ~= 0?
        fixQ,               ...   % Condition on Q ~= 10^{-5}?
        Nafter,             ...   % Number of e-folds between phiexit and phiend
        seed,               ...   % Seed for the random number generator
        n_recycle           ...   %
    )
    
    if nargin < 1,  outfile = '';           end
    if nargin < 2,  n_iter = 1e3;           end
    if nargin < 3,  mv = 0.007;             end
    if nargin < 4,  mh = 1;                 end
    if nargin < 5,  kmax = 30;              end
    if nargin < 6,  gamma = 0;              end
    if nargin < 7,  measure = 'B';          end
    if nargin < 8,  n_tunnel_max = 3;       end
    if nargin < 9,  lambdascreen = true;    end
    if nargin < 10, rho_Lambda_thres = 1e-7;end
    if nargin < 11, fixQ = false;           end
    if nargin < 12, Nafter = 55;            end
    if nargin < 13 || seed == -1, seed = randi(1e5); end
    if nargin < 14, n_recycle = 4;          end
    
%     seed = 94335;

    %% Create output file
    
%     fid = fopen(outfile,'a');
%     if fid == -1
%         fid = fopen(outfile,'w'); fclose(fid);
%         fid = fopen(outfile,'a');
%     end
%     fprintf(fid,'%E,%.4G,%.4G,%d,%.2f,%s,%d,%d,%.4G,%d,%.2f\r\n',...
%         n_iter,mv,mh,kmax,gamma,measure,n_tunnel_max,...
%         logical(lambdascreen),rho_Lambda_thres,logical(fixQ),Nafter);
%     fclose(fid);
    
    disp(['outfile: ' outfile]);
    
    %% Initialize and run eternal inflation simulator
    
    eis = EternalInflationSimulator(...
        'seed',               seed,...
        'n_iter',             n_iter,...
        'mv',                 mv,...
        'mh',                 mh,...
        'kmax',               kmax,...
        'gamma',              gamma,...
        'Nafter',             Nafter,...
        'lambdascreenmode',   logical(lambdascreen),...
        'rho_Lambda_thres',   rho_Lambda_thres,...
        'fixQ',               logical(fixQ),...
        'measure',            measure,...
        'n_tunnel_max',       n_tunnel_max,...
        'outfile',            outfile,...
        'n_recycle',          n_recycle );
    
    eis.main()
    
    
    
    
    