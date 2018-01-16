n_iter = 5e3;

mv_range = [0.007];
mh_range = [1];

var_ranges = {mv_range,mh_range};
n_test = prod(cellfun(@length,var_ranges));
n_vars = length(fieldnames(EternalInflationSimulator.results_template));

results = sparse(n_test*n_iter,n_vars);

% %% Load params from a log file
% 
% fnames = dir('../data');
% fnames = fnames(3:end);
% fnames = {fnames([~fnames.isdir]).name};
% 
% numdatasets = length(fnames);
% 
% for i = 1:numdatasets
%     
%     fname = fnames{i};
%     
% %%

eis = EternalInflationSimulator(...
    'n_iter',n_iter);

objects = cell(n_test,1);

i_test = 0;
for mv = mv_range
    for mh = mh_range
        
        i_test = i_test + 1;
%         i_iter = i_iter+n_iter;
        
        eis = EternalInflationSimulator(...
            'n_iter',n_iter,'mv',mv,'mh',mh,'fixLambda',false);
        eis.main()
        
        objects{i_test} = eis;
        
%         results(n_iter*(i_test-1)+(1:n_iter),:) = eis.results;
        
    end
end

save(['eternal_sim_outputs_' date()],'objects')
