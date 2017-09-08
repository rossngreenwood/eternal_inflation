eis = EternalInflationSimulator();

n_iter = 1e3;

mv_range = 1:10;
mh_range = 1:10;

var_ranges = {mv_range,mh_range};
n_test = sum(cellfun(@length,var_ranges));
n_vars = length(fieldnames(EternalInflationSimulator.results_template));

outputs = single(zeros(n_test*n_iter,n_vars));

i_test = 1-n_iter;
for mv = mv_range
    for mh = mh_range
        
        i_test = i_test+n_iter;
        
        eis.set_parameters('mv',mv,'mh',mh)
        eis.main()
        
        outputs(i_test-n_iter:i_test,:) = eis.results;
        
    end
end

save eternal_sim_outputs outputs -append
