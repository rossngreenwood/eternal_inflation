classdef EternalInflationSimulator < handle
% Running EternalInflationSimulator.main() performs a set of inflation
% simulations and records data relevant to the presence of eteranl
% inflation. At each iteration, the code generates a potential as a
% Gaussian random field and initializes the inflaton field according to
% some specified measure.

%% Properties, constructor, and set/get methods

properties
    
    parameters = struct(...
        'n_iter',             1e4,...   % Number of iterations
        'mv',                 1,...     % Initial mass scale of the potential
        'mh',                 5,...     % Mass scale of the inflaton field
        'kmax',               50,...    % Largest wavenumber for GRF
        'gamma',              0,...     % Frequency dependence of GRF
        'Nafter',             55,...    % Number of e-folds between phiexit and phiend
        'lambdascreenmode',   true,...  % Throw out cases where rho_Lambda < 0?
        'fixLambda',          false,... % Condition on rho_Lambda ~= 0?
        'fixQ',               false,... % Condition on Q ~= 10^{-5}?
        'measure',            'B',...   % Measure on initial conditions
        'n_tunnel_max',       1);       % Max number of tunneling events to simulate
    
    results     % Table of observables
    
end

properties (Constant)
    
    results_template = struct(...
            'status',        SRStatus.null,...  % Simulation status
            ...
            'Ntotal',        nan('single'),...  % Observables
            'Q',             nan('single'),...
            'r',             nan('single'),...
            'n_s',           nan('single'),...
            'alpha',         nan('single'),...
            'n_t',           nan('single'),...
            'dlgrho',        nan('single'),...
            'lgOk',          nan('single'),...
            'rho_Lambda',    nan('single'),...
            ...
            'phitunnel',     nan('single'),...  % False-vacuum eternal
            'tunneling_rate',nan('single'),...
            ...
            'NStochastic',   nan('single'),...  % Stochastic eternal
            'NSinceStoch',   nan('single'),...
            'numStochEpochs',int8(-1),...
            ...
            'numTopolEpochs',int8(-1),...       % Topological eternal
            'mv',            1);
    
end

methods
    
    function obj = EternalInflationSimulator(params,varargin)
        % Constructor
        % Inputs (varargin)
        % 1) params     A structure containing some or all of the fields in
        %               EternalInflationSimulator.parameters
        % 2) {field1,val1,field2,val2,...}  Field and value pairs
        
        if nargin < 1, params = struct(); end
        
        obj.set_parameters(params,varargin{:});
        
    end
    
    function set_parameters(obj,varargin)
        % Set simulation parameters
        % Inputs (varargin)
        % 1) params     A structure containing some or all of the fields in
        %               EternalInflationSimulator.parameters
        % 2) {field1,val1,field2,val2,...}  Field and value pairs
        
        if nargin == 2
            val = varargin{:};
        else
            val = struct(varargin{:});
        end
        
        for fn = fieldnames(val).'
            switch fn{:}
                
                case 'n_iter'
                    
                    if isnumeric(val.n_iter) && isscalar(val.n_iter) && ...
                            isfinite(val.n_iter) && val.n_iter > 0 && ...
                            mod(val.n_iter,1) == 0 && isreal(val.n_iter)
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.n_iter must be a real, finite, positive integer.');
                    end
                    
                case 'mv'
                    
                    if isnumeric(val.mh) && isscalar(val.mh) && ...
                            isfinite(val.mh) && val.mh > 0 && isreal(val.mv)
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.mv must be a real, finite, positive number.');
                    end
                    
                case 'mh'
                    
                    if isnumeric(val.mh) && isscalar(val.mh) && ...
                            isfinite(val.mh) && val.mh > 0 && isreal(val.mh)
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.mh must be a real, finite, positive number.');
                    end
                    
                case 'kmax'
                    if isnumeric(val.kmax) && isscalar(val.kmax) && ...
                            isfinite(val.kmax) && isreal(val.kmax)
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.kmax must be a real, finite number.');
                    end
                    
                case 'gamma'
                    if isnumeric(val.gamma) && isscalar(val.gamma) && ...
                            isfinite(val.gamma) && isreal(val.gamma)
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.gamma must be a real, finite number.');
                    end
                    
                case 'Nafter'
                    if isnumeric(val.Nafter) && isscalar(val.Nafter) && ...
                            isfinite(val.Nafter) && val.Nafter > 0 && ...
                            isreal(val.Nafter)e
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.Nafter must be a finite, positive number.');
                    end
                    
                case 'lambdascreenmode'
                    if islogical(val.lambdascreenmode) && isscalar(val.lambdascreenmode)
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.lambdascreenmode must be a logical.');
                    end
                    
                case 'fixLambda'
                    if islogical(val.fixLambda) && isscalar(val.fixLambda)
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.fixLambda must be a logical.');
                    end
                    
                case 'fixQ'
                    if islogical(val.fixQ) && isscalar(val.fixQ)
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.fixQ must be a logical.');
                    end
                    
                case 'measure'
                    if ismember(upper(val.measure),{'A','B'})
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.measure must be ''A'' or ''B''');
                    end
                    
                case 'n_tunnel_max'
                    if isnumeric(val.n_tunnel_max) && isscalar(val.n_tunnel_max) && ...
                            isfinite(val.n_tunnel_max) && val.n_tunnel_max > 0 && ...
                            mod(val.n_tunnel_max,1) == 0 && isreal(val.n_tunnel_max)
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.n_tunnel_max must be a real, finite, positive integer.');
                    end
                    
            end
        end
    end
    
end

%% Main method

methods (Access = public)
    
    function main(obj)
        % Collect data from inflation simulations pertaining to the
        % manifestation of false-vacuum, stochastic, and topological
        % eternal inflation, as well as CMB observables.
        
        p = obj.parameters;
        
        % Use sparse data matrix to reduce memory demands
        obj.results = sparse(p.n_iter,numel(fieldnames(obj.results_template)));
        
        % Run n_iter iterations of the inflation simulation
        for i_iter = 1:p.n_iter
            
            for goto = NaN % on break, goto record data
                
                data = obj.results_template;
                
                %% Set initial conditions
                
                [~,f{1:4}] = obj.gaussian_random_field_1D(p.kmax,p.gamma);
                [data_init,phistart,flagStartAtMaximum] = obj.set_initial_conditions(f);
                data = structcat(data,data_init);
                
                %% Simulate slow roll inflation
                
                [data_sr,pot,phi] = obj.simulate_slowroll(...
                    f,phistart,flagStartAtMaximum);
                
                if isnan(data_sr.Ntotal)
                    break % No inflation; skip to record data
                end
                
                %% Check for false-vacuum eternal inflation
                
                for n_tunnel = 1:p.n_tunnel_max
                    
                    % Look for an instanton tunneling solution
                    [data_fv] = obj.check_false_vacuum_eternal(pot,phi);
                    data = structcat(data,data_fv);
                    
                    if isnan(data_fv.phitunnel)
                        break % No tunneling; move on
                    end
                    
                    % If tunneling occurs, simulate slow roll from the new
                    % starting point on other side of potential barrier
                    [data_sr,pot,phi] = obj.simulate_slowroll(f,data.phitunnel);
                    
                end
                
                % Keep only the last bout of inflation
                data = structcat(data,data_sr);
                
                if data.Ntotal < obj.parameters.Nafter
                    break % Not enough inflation; skip to record data
                end
                
                %% Check for stochastic and topological eternal inflation
                
                data_stoch = obj.check_stochastic_eternal(pot,phi,data_sr);
                data_topol = obj.check_topological_eternal(pot,phi);
                data = structcat(data,data_stoch,data_topol);
                
                %% Compute observables
                
                Nbefore = data_sr.Ntotal - p.Nafter; % e-folds before crossing
                observables = obj.compute_observables(pot,phi,Nbefore);
                data = structcat(data,observables);
                
            end % for goto
            
            %% Record results
            
            data = struct2cell(data);
            data = cellfun(@double,data,'Un',0);
            data = horzcat(data{:});
            data(isnan(data)) = 0;
            
            obj.results(i_iter,:) = data;
            
        end % for i_iter
        
    end
    
end

%% Slowroll

methods (Access = protected)
    
    function [datastruct,phistart,startAtMaximum] = set_initial_conditions(obj,f)
        % Simulate slow roll inflation and collect data
        %
        % Inputs
        %   f   A cell array of function handles
        %           {V(phi), V'(phi), V''(phi), V'''(phi)}
        %
        % Outputs
        %   datastruct  A structure containing observables and diagnostics
        
        p = obj.parameters;
        
        datastruct = obj.results_template;
        datastruct.mv = p.mv;
        
        [V,Vp,Vpp] = build_potential(f,p.mv,p.mh);
        
        %% Initialize field value according to Measure A/B
        
        switch upper(p.measure)
            case 'A'
                % Start at maximum (first slowroll region below max)
                if V(0) < 0
                    datastruct.status = 5;
                    return
                end
                [phistart,status] = obj.find_phistart_A(0,V,Vp,Vpp,p.mh);
                if status == 4 && phistart ~= 0
                    phistart = phistart + 1i; % Start at maximum
                end
            case 'B'
                phistart = 0;
        end
        
    end
    
    function [datastruct,potential,phi] = simulate_slowroll(obj,f,phistart)
        % Simulate slow roll inflation and collect data
        %
        % Inputs
        %   f   A cell array of function handles
        %           {V(phi), V'(phi), V''(phi), V'''(phi)}
        %
        % Outputs
        %   datastruct  A structure containing observables and diagnostics
        
        startAtMaximum = ~isreal(phistart);
        phistart = real(phistart);
        
        p = obj.parameters;
        
        datastruct = obj.results_template;
        datastruct.mv = p.mv;
        
        phi = nan(1,5);
        phi(2) = phistart;
        
        mh = p.mh;
        
        %% Build potential
        
        [V,Vp,Vpp,Vppp] = build_potential(f,p.mv,p.mh);
        potential = {V,Vp,Vpp,Vppp};
        
        %% Simulate slow roll
        
        % Slow roll conditions
        slowroll = @(x) chop((Vp(x)./V(x)).^2/2,1) & chop(abs(Vpp(x)./V(x)),1);
        
        shift_cases = {'Unchanged'};
        if p.fixLambda, shift_cases = [shift_cases {'Lambda'}]; end
        if p.fixQ,      shift_cases = [shift_cases {'Q'}];      end
        
        more_inflation = true;
        while more_inflation
            
            for shifted_potential = shift_cases
                
                %% Modify potential to match observables
                
                switch shifted_potential{:}
                    
                    case 'Lambda'
                        
                        % If present-day cosmological constant is above some
                        % threshold, subtract it from the potential and retry
                        thres = 1e-20;
                        V_min = V(phistop);
                        if abs(V_min) > thres
                            A = prefactor{1}; grf = f{1};
                            V = @(phi) A*grf(phi/mh) - V_min;
                        end
                        
                    case 'Q'
                        
                        % Rescale mv so that Q matches observation
                        Q_target = 2e-5;
                        Q = sqrt(V(phiexit)/(150*(Vp(phiexit)./V(phiexit)).^2/2))/pi;
                        datastruct.mv = obj.parameters.mv*sqrt(Q_target/Q);
                        
                        [V,Vp,Vpp] = build_potential(f,sqrt(Q_target/Q)*p.mv,p.mh);
                        
                end
                
                %% Identify phiend, phiexit, and phistop
                
                % Find the end of slow roll inflation
                [phiend,datastruct.status] = obj.find_phiend(phistart,V,Vp,Vpp,p.mh);
                phi(4) = phiend;
                if phiend == phistart, break, end
                if datastruct.status == 4, return, end
                
                % Throw out models that settle into a negative minimum
                phistop = obj.find_phistop(phiend,V,Vp,Vpp,mh);
                phi(5) = phistop;
                if p.lambdascreenmode
                    if V(phistop) < 0
                        datastruct.status = datastruct.status + 5;
                        return
                    end
                end
                
                % Compute total number of slow roll e-folds
                dlna_dphi = @(phi) (-V(phi)./Vp(phi));
                switch upper(p.measure)
                    case 'A'
                        if startAtMaximum
                            datastruct.Ntotal = Inf;
                        else
                            datastruct.Ntotal = integral(dlna_dphi,phistart,phiend);
                        end
                    case 'B'
                        datastruct.Ntotal = integral(dlna_dphi,phistart,phiend);
                end
                if datastruct.Ntotal < p.Nafter, return, end
                
                % Find the value of the field 55 e-folds
                % before the end of slow roll inflation
                dphi_dlna = @(N,phi) (-Vp(phi)./V(phi)).';
                Nbefore = datastruct.Ntotal - p.Nafter;
                [~,phiexit] = ode45(dphi_dlna,[0 Nbefore Nbefore+1],phistart);
                phiexit = phiexit(2);
                phi(3) = phiexit;
                
                % Find the nearest minimum after inflation ends
                if ~p.lambdascreenmode
                    phistop = obj.find_phistop(phiend,V,Vp,Vpp,mh);
                end
                phi(5) = phistop;
                
            end
            
            phimax = fminsearch(@(x) -V(x),phistart);
            phi(1) = phimax;
            
            if phistart == phiend, return, end
            
            break
            
            %% Look beyond phiend for more slow roll inflation
            
            % Nearest maximum
            phimax = fminsearch(@(x) -V(x),phistart);
            
            philower = phistop + sgn(Vp(phistart))*0.01*mh;
            Nlower = integral(dlna_dphi,phiend,philower);
            
            lna_steps = 0:floor(Nlower);
            [~,phi_N] = ode45(dphi_dlna,lna_steps,phistart);    % Phi values separated by 1 e-fold of inflation
            if phimax > phistart, phi_N = fliplr(phi_N); end   	% Order [phimax; ... phistart; ... phiend]
            
            phi_sr = phi_N(slowroll(V(phi_N)));
            more_inflation = ~isempty(phi_sr);
            
            look_for_future_inflation();
            sample_phidot();
            do_full_simulation();
            
        end
        
        % Move phistart to just shy of the local maximum
        if strcmpi(p.measure,'A')
            delta = 1e-8*p.mh;
            if Vp(phistart) > 0, delta = -delta; end
            phistart = fminsearch(@(x) -V(x),phistart) + delta; 
        end
        
        phi = [phimax,phistart,phiexit,phiend,phistop];
        
    end
    
end

%% Checks for eternal inflation

methods (Access = protected)
    
    function [datastruct] = check_false_vacuum_eternal(obj,potential,phi)
        % Handle false vacuum tunneling.
        
        p = obj.parameters;
        
        [V,Vp,Vpp] = deal(potential{1:3});   % Potential function
        phistop = phi(end); % Location of minimum of basin
        
        %% Locate nearest minima
        
        % Find minima around phistop
        phispace = phistop + linspace(-10*p.mh,10*p.mh,1001);
        [~,approx_mins] = findpeaks(-V(phispace));
        [~,imin] = min(abs(phispace(approx_mins)-phistop));
        
        % Find nearest minima
        near_minima = nan(1,2);
        if 1 < imin
            near_minima(1) = phispace(approx_mins(imin-1));
        end
        if imin < length(approx_mins)
            near_minima(2) = phispace(approx_mins(imin+1));
        end
        near_minima(V(near_minima) >= V(phistop)) = NaN;
        
        %% Check for tunneling to nearby minima
        
        new_phistart   = nan(size(near_minima));
        tunneling_rate = zeros(size(near_minima));
        for lr = 1:2
            
            %% Try less intensive methods for computing tunneling rate
            %  Thin wall approximation valid?
            
            do_rigorous_calculation = true;
            
            if ~do_rigorous_calculation, continue, end
            
            %% If approximate methods fail, do full instanton calculation
            
            if isnan(near_minima(lr)), continue, end
            
            fvi = FalseVacuumInstanton(...
                'V',            V,...
                'dV',           Vp,...
                'd2V',          Vpp,...
                'M_Pl',         sqrt(8*pi),...
                'phi_metaMin',  phistop, ...
                'phi_absMin',   near_minima(lr) );
            
            %% Solve for instanton profile
            
            % Defaults
            xguess           = [];
            xtol             = 1e-4;
            phitol           = 1e-4;
            thinCutoff       = 0.01;
            npoints          = 500;
            rmin             = 1e-4;
            rmax             = 1e4;         % Units?
            max_interior_pts = [];
            
            [R,Y,~,n_interior_pts] = fvi.find_profile(xguess,xtol,phitol,...
                thinCutoff,npoints,rmin,rmax,max_interior_pts);
            
            %% Compute tunneling rate
            
            % Solution exterior to bubble
            R_ext = R(n_interior_pts+1:end);
            Y_ext = Y(:,n_interior_pts+1:end);
            
            new_phistart(lr) = Y(1,1); % field value at center of bubble
            tunneling_rate(lr) = fvi.find_tunneling_rate(R,Y);
            
        end
        
        % Choose bubble with the larger tunneling rate
        [datastruct.tunneling_rate,imax] = max(tunneling_rate);
        
        rate_thres = 9/4/pi * H^4; % Units? Justification?
        if tunneling_rate >= rate_thres
            % Find the new value of phi after tunneling
            datastruct.phitunnel = new_phistart(imax);
        else
            datastruct.phitunnel = NaN;
        end
        
    end
    
    function [datastruct] = check_stochastic_eternal(obj,potential,phi,datastruct_sr)
        
        p = obj.parameters;
        
        datastruct = struct();
        
        if isnumeric(phi), phi = num2cell(phi); end
        
        [phimax,phistart,phiexit,phiend] = deal(phi{1:4});
        [V,Vp,Vpp] = deal(potential{1:3});
        
        dlna_dphi = @(phi) (-V(phi)./Vp(phi));
        dphi_dlna = @(N,phi) (-Vp(phi)./V(phi)).';
        slowroll = @(x) chop((Vp(x)./V(x)).^2/2,1) & chop(abs(Vpp(x)./V(x)),1);
        
        %% Stochastic Eternal inflation
        
        stochasticEIC = @(x) V(x).^(3/2) > 2*pi*sqrt(3)*(0.61)*abs(Vp(x)); % stochasticEIC > 0 -> Eternal
        
        % Find the value of the field 55 e-folds
        % before the end of slow roll inflation
        
        lna_steps = 0:floor(datastruct_sr.Ntotal);
        if length(lna_steps) == 1, return, end
        [~,phi_N] = ode45(dphi_dlna,lna_steps,phistart);                 % Phi values separated by 1 e-fold of inflation
        phi_N = sort([phi_N; phiend; linspace(phistart,phimax,20).']).'; % Add search points between phistart and phimax
        if phimax > phistart, phi_N = fliplr(phi_N); end                 % Order [phimax; ... phistart; ... phiend]
        phi_N(isnan(phi_N)) = [];
        
        % Find where stochastic eternal inflation starts/ends
        % Should always be eternal at phimax
        intv_ind = find(diff(stochasticEIC(phi_N)));
        off2on_sei   = feval(@(x) x(intv_ind) > 0, diff(stochasticEIC(phi_N)));
        phibreak_sei = arrayfun(@(i) fzero(stochasticEIC,phi_N(i:i+1)),intv_ind);
        
        disp(num2str(length(phibreak_sei)));
        
        disp([num2str(V(phistart)^(3/2)) ' ' num2str(2*pi*sqrt(3)*(0.61)*abs(Vp(phistart)))]);
        
        if isempty(phibreak_sei)
            phibreak_sei = phimax;
            off2on_sei = 0;
        end
        
        % Second derivative check for stochastic eternal inflation near the
        % maximum. Assume potential is locally mirror symmetric about maximum,
        % or we condition on the field falling only toward phistart.
        phidev = phimax-phibreak_sei(1); Hmax = sqrt(V(phimax)/3);
        if abs(erf(abs(phidev/(Hmax/2/pi))/sqrt(2))) > exp(-3)
            datastruct.numStochEpochs = uint8(1 + nnz(off2on_sei));
        else
            stochasticEIC = @(x) stochasticEIC(x) & ...
                ~(min(phimax,phibreak_sei(1)) <= x && x <= max(phimax,phibreak_sei(1)));
            phibreak_sei(1) = []; off2on_sei(1) = [];
            datastruct.numStochEpochs = uint8(nnz(off2on_sei));
        end
        
        % Compute # of e-folds after past SEI breakdown and before phiexit
        if stochasticEIC(phiexit)
            datastruct.NSinceStoch = 0; % Eternal at phiexit
        elseif any((phibreak_sei(~off2on_sei)-phiexit)*Vp(phiexit) > 0)
            phib = phibreak_sei((phibreak_sei-phiexit)*Vp(phiexit) > 0 & ~off2on_sei);
            [~,imin] = min(abs(phib-phiexit));
            datastruct.NSinceStoch = integral(@(x) dlna_dphi(x).*slowroll(x),...
                phib(imin),phiexit);
        else
            datastruct.NSinceStoch = NaN; % Never eternal before phiexit
        end
        
        % Compute # of e-folds eternal after phistart and before phiend
        if strcmpi(p.measure,'a') && stochasticEIC(phistart)
            datastruct.NStochastic = Inf;
        else
            datastruct.NStochastic = integral(@(x) ...
                dlna_dphi(x).*stochasticEIC(x).*slowroll(x),phistart,phiend);
        end
        
    end
    
    function [datastruct] = check_topological_eternal(obj,potential,phi)
        
        p = obj.parameters;
        
        datastruct = obj.results_template;
        
        if isnumeric(phi), phi = num2cell(phi); end
        
        [V,Vp,Vpp] = deal(potential{1:3});
        [phimax,phistart] = deal(phi{1:2});
        
        %% Topological Eternal Inflation
        
        % Determine whether quantum fluctuations could result in at least
        % one Hubble volume descending toward a different minimum of V(phi)
        if strcmpi(p.measure,'a')
            checkTopological = true;
        else
            H = sqrt(V(phistart)/3); dphi = H/2/pi; Dphi = -Vp(phistart)/(3*H^2);
            checkTopological = 1/2*erfc(abs(phistart+Dphi-phimax)/(sqrt(2)*dphi)) > exp(-3);
        end
        
        if checkTopological
            
            phimax = phi{1};
            
            % Find value of phi at the domain wall boundary
            phiedge_eps = fzero(@(phi) Vp(phi)./V(phi)/sqrt(2) - 1,phimax);
            phiedge_eta = fzero(@(phi) abs(Vpp(phi)./V(phi)) - 1,phimax);
            [~,icloser] = min(abs(phimax-[phiedge_eps,phiedge_eta]));
            phiedge = feval(@(x) x(icloser),[phiedge_eps,phiedge_eta]);
            
            % Find value of phi that will descend to phistart in time < t_H
            phistar = NaN; q = 1;
            while isnan(phistar)
                try
                    phistar = fzero(@(phi) phi + Vp(phi)./V(phi) - phiedge, [phimax q*phiedge]);
                catch ME
                    if strcmp(ME.identifier,'MATLAB:fzero:ValuesAtEndPtsSameSign')
                        q = q + 1;
                        continue
                    end
                    rethrow(ME);
                end
            end
            
            % If the expansion of the inner portion of the domain wall
            % where (phimax < phi < phistar) replaces loss of the outer
            % wall where (phistar < phi < edge), then inflation is eternal
            if abs(phistar-phimax) > abs(phiedge-phimax)*exp(-1)
                datastruct.numTopolEpochs = true;
            end
            
        end
        
    end
    
end

methods
    
end

%% Subroutines

methods (Static, Access = protected)
    
    function look_for_future_inflation(phistart,V,Vp,Vpp,phiscale)
        
    end
    
    function [phiend,status] = find_phiend(phistart,V,Vp,Vpp,phiscale)
        % Find the value of phi at the end of slow roll inflation
        
        status = SRStatus.compute_status(V,Vp,Vpp,phistart);
        if ~status.isok(), phiend = phistart; return, end % DOA
        
        sgn    = sign(Vp(phistart));
        phimin = phiscale^2 * abs(Vp(phistart)./V(phistart));
        dphi   = -0.001*sgn*max(0.01*phiscale,min(phimin,phiscale));            %%% Tunable
        
        % Take steps until passed breakdown of slow roll
        step = 10^(1/16);                                                       %%% Tunable
        while status.isok()
            phiend = phistart + dphi;
            status = SRStatus.compute_status(V,Vp,Vpp,phiend);
            if sgn*Vp(phiend) < 0, status = SRStatus.localmin; end
            dphi = dphi*step;
        end
        
        % Choose function to zero based on status
        switch status
            case SRStatus.rho, fun = V;                            % V = 0
            case SRStatus.eps, fun = @(x) 1 - (Vp(x)./V(x)).^2/2;  % eps = 1
            case SRStatus.eta, fun = @(x) 1 - abs(Vpp(x)./V(x));   % eta = 1
            case SRStatus.localmin, fun = Vp;                           % Vp = 0
        end
        
        % Find zero of function to obtain true phiend
        phiend = fzero(fun,[phiend phistart]);
        
    end
    
    function [phistop] = find_phistop(phiend,V,Vp,~,phiscale)
        % Find the value of phi at the next local minimum
        
        sgn    = sign(Vp(phiend));
        phimin = phiscale^2 * abs(Vp(phiend)./V(phiend));
        dphi   = -0.001*sgn*max(0.01*phiscale,min(phimin,phiscale));            %%% Tunable
        
        % Take steps until passed local min
        phistop = phiend;
        step = 10^(1/16);                                                       %%% Tunable
        while sgn*Vp(phistop) > 0
            phistop = phiend + dphi;
            dphi = dphi*step;
        end
        
        % Find zero of V'(phi) to obtain true min
        phistop = fzero(Vp,[phistop phiend]);
        
    end
    
    function [phistart,status] = find_phistart_A(phi0,V,Vp,Vpp,phiscale)
        % Find the value of phi at the end of slow roll inflation
        
        % Screen negative starting points?
        
        orig_status = SRStatus.compute_status(V,Vp,Vpp,phi0);
        status = orig_status;
        
        sgn    = sign(Vp(phi0));
        phimin = phiscale^2 * abs(Vp(phi0)./V(phi0));
        dphi   = -0.001*sgn*max(0.01*phiscale,min(phimin,phiscale));            %%% Tunable
        
        sr_at_0 = status.isok(); % Is SRA valid at phistart?
        if sr_at_0
            search_direction = -1; % If inflating, look up
        else
            search_direction = 1;  % If not inflating, look down
        end
        
        % Take steps until passed breakdown of slow roll
        step = 10^(1/16);                                                       %%% Tunable
        while (sr_at_0 && status.isok()) || (~sr_at_0 && status.isok())
            phistart = phi0 + search_direction*dphi;
            status = SRStatus.compute_status(V,Vp,Vpp,phistart);
            if sgn*Vp(phistart) < 0
                status = SRStatus.localmin;
                if sr_at_0
                    phistart = phi0;
                    return % No SR before min 
                end
            end
            dphi = dphi*step;
        end
        
        if ~sr_at_0, status = orig_status; end
        
        % Choose function to zero based on status
        switch status
            case SRStatus.rho, fun = V;                            % V = 0
            case SRStatus.eps, fun = @(x) 1 - (Vp(x)./V(x)).^2/2;  % eps = 1
            case SRStatus.eta, fun = @(x) 1 - abs(Vpp(x)./V(x));   % eta = 1
            case SRStatus.localmin, fun = Vp;                      % Vp = 0
        end
        
        % Find zero of function to obtain true phistart
        phistart = fzero(fun,[phistart phi0]);
        
    end
    
    function [status] = compute_status(V,Vp,Vpp,phi)
        if     V(phi) < 0,                status = 1;     % \rho_Lambda < 0
        elseif (Vp(phi)/V(phi)).^2/2 > 1, status = 2;     % \eps > 1
        elseif abs(Vpp(phi)/V(phi)) > 1,  status = 3;     % \eta > 1
        else                              status = 0; end % OK
    end
    
end

methods (Static)
    
    function [observables] = compute_observables(potential,phi,Nbefore)
        % Inputs
        %   potential   Cell array of function handles
        %                   {V, V', V'', V'''}
        %   phi         Monotonic array of field values
        %                   [phistart, phiexit, phiend, phistop]
        %   Nbefore     Number of e-foldings between phistart and phiexit
        %
        % Outputs
        %   observables Structure of observables
        
        if isnumeric(phi), phi = num2cell(phi); end
        
        [phistart,phiexit,phiend,phistop] = deal(phi{1:4});
        [V,Vp,Vpp,Vppp] = deal(potential{1:4});
        
        V_exit = V(phiexit);
        
        % Slow roll parameters at horizon exit scale
        eps = (Vp(phiexit)./V_exit).^2/2;
        eta = Vpp(phiexit)./V_exit;
        xi2 = Vp(phiexit).*Vppp(phiexit)./V_exit^2;
        
        observables.Q           = sqrt(V_exit/(150*eps))/pi;
        observables.r           = 16*eps;
        observables.n_s         = 1-6*eps+2*eta;
        observables.alpha       = 16*eps*eta - 24*eps^2 - 2*xi2;
        observables.n_t         = -2*eps;
        observables.dlgrho      = log10(V_exit/V(phiend));
        observables.lgOk        = log10(V(phistart)/V_exit) - Nbefore*2/log(10);
        observables.rho_Lambda  = V(phistop);
        
    end
    
end

%% Utilities

methods (Static)
    
    function plot_histograms(observables,savename)
        
        if nargin < 2, savename = ''; end
        save_fig = [savename '_%s'];
        
        save(sprintf(save_fig,'observables'),'observables');
        
        style = hgexport('factorystyle');
        style.Height = 6;
        style.Width  = 7;
        style.Bounds = 'tight';
        style.FontSizeMin = 20;
        style.LineWidthMin = 1;
        
        hist_color = [142 219 255]/255.0;
        
        % Power
        fieldname = 'Q';
        figure, hist(log10([observables.(fieldname)]),50,hist_color,'Linewidth',0)
        title('Power: Q'), set(gca,'XTick',-6:4)
        set(gca,'XTickLabel',arrayfun(@(x) sprintf('10^{%1i}',x),-6:4,'Un',0))
        h = findobj(gca,'Type','patch'); h.FaceColor = hist_color; boldify
        if ~isempty(savename)
            hgexport(gcf,'-clipboard',style,'applystyle',true); drawnow
            savefig(sprintf(save_fig,fieldname));
            print(sprintf(save_fig,fieldname),'-dpng','-r0');
        end
        close(gcf)
        
        % n_s - Scalar spectral index
        fieldname = 'n_s';
        figure, hist([observables.(fieldname)],-0.01:0.05:2,hist_color)
        title('Scalar Spectral Index: n_s')
        set(gca,'XLim',[0 2]);
        h = findobj(gca,'Type','patch');
        h.FaceColor = hist_color; boldify
        if ~isempty(savename)
            hgexport(gcf,'-clipboard',style,'applystyle',true); drawnow
            savefig(sprintf(save_fig,fieldname));
            set(gcf,'PaperPositionMode','auto');
            print(sprintf(save_fig,fieldname),'-dpng','-r0');
        end
        close(gcf)
        
        % alpha - Running of the spectral index
        fieldname = 'alpha';
        figure, hist([observables.(fieldname)],-0.51:0.02:0.51,hist_color)
        title('Running: \alpha')
        set(gca,'XLim',[-0.5,0.5]);
        h = findobj(gca,'Type','patch');
        h.FaceColor = hist_color; boldify
        if ~isempty(savename)
            hgexport(gcf,'-clipboard',style,'applystyle',true); drawnow
            savefig(sprintf(save_fig,fieldname));
            set(gcf,'PaperPositionMode','auto');
            print(sprintf(save_fig,fieldname),'-dpng','-r0');
        end
        close(gcf)
        
        % n_t - Tensor spectral index
        fieldname = 'n_t';
        figure, hist([observables.(fieldname)],-0.051:0.001:0.01,hist_color)
        title('Tensor Spectral Index: n_t')
        set(gca,'XLim',[-0.05,0]);
        h = findobj(gca,'Type','patch');
        h.FaceColor = hist_color; boldify
        if ~isempty(savename)
            hgexport(gcf,'-clipboard',style,'applystyle',true); drawnow
            savefig(sprintf(save_fig,fieldname));
            set(gcf,'PaperPositionMode','auto');
            print(sprintf(save_fig,fieldname),'-dpng','-r0');
        end
        close(gcf)
        
        % log |Omega_tot - 1|
        fieldname = 'lgOk';
        figure, hist(log10(abs([observables.(fieldname)])),50,hist_color)
        title('$$\log_{10} \vert \Omega_{tot} - 1 \vert $$','Interpreter','LaTex'), set(gca,'XTick',1:5)
        set(gca,'XTickLabel',arrayfun(@(x) sprintf('10^{%1i}',x),1:5,'Un',0))
        set(gca,'XLim',[0,5]);
        h = findobj(gca,'Type','patch');
        h.FaceColor = hist_color; boldify
        if ~isempty(savename)
            hgexport(gcf,'-clipboard',style,'applystyle',true); drawnow
            savefig(sprintf(save_fig,fieldname));
            set(gcf,'PaperPositionMode','auto');
            print(sprintf(save_fig,fieldname),'-dpng','-r0');
        end
        close(gcf)
        
        % Number of e-foldings since eternal inflation
        fieldname = 'NSinceStoch';
        val = log10([observables.(fieldname)]);
        figure, hist(val(isreal(val)),50,hist_color)
        title('Number of e-foldings since eternal');
        xlim = get(gca,'XLim');
        set(gca,'XTick',ceil(xlim(1)):floor(xlim(2)))
        set(gca,'XTickLabel',arrayfun(@(x) sprintf('10^{%1i}',x),...
            ceil(xlim(1)):floor(xlim(2)),'Un',0))
        h = findobj(gca,'Type','patch'); h.FaceColor = hist_color; boldify
        if ~isempty(savename)
            hgexport(gcf,'-clipboard',style,'applystyle',true); drawnow
            savefig(sprintf(save_fig,fieldname));
            set(gcf,'PaperPositionMode','auto');
            print(sprintf(save_fig,fieldname),'-dpng','-r0');
        end
        close(gcf)
        
        % Number of e-foldings of eternal inflation
        fieldname = 'NStochastic';
        q = log10([observables.(fieldname)]);
        figure, hist(q(~isinf(q)),50,hist_color)
        title('Number of e-foldings eternal');
        xlim = get(gca,'XLim');
        set(gca,'XTick',ceil(xlim(1)):floor(xlim(2)))
        set(gca,'XTickLabel',arrayfun(@(x) sprintf('10^{%1i}',x),...
            ceil(xlim(1)):floor(xlim(2)),'Un',0))
        h = findobj(gca,'Type','patch'); h.FaceColor = hist_color; boldify
        if ~isempty(savename)
            hgexport(gcf,'-clipboard',style,'applystyle',true); drawnow
            savefig(sprintf(save_fig,fieldname));
            set(gcf,'PaperPositionMode','auto');
            print(sprintf(save_fig,fieldname),'-dpng','-r0');
        end
        close(gcf)
        
        % Number of eternal inflation epochs
        fieldname = 'numStochEpochs';
        figure, hist([observables.(fieldname)],10,hist_color)
        title('Number of eternal inflation epochs')
        h = findobj(gca,'Type','patch');
        h.FaceColor = hist_color; boldify
        if ~isempty(savename)
            hgexport(gcf,'-clipboard',style,'applystyle',true); drawnow
            savefig(sprintf(save_fig,fieldname));
            set(gcf,'PaperPositionMode','auto');
            print(sprintf(save_fig,fieldname),'-dpng','-r0');
        end
        close(gcf)
        
    end
    
    function [a,varargout] = gaussian_random_field_1D(n,gamma,a)
        % [a,V,Vp,Vpp,...] = gaussian_random_field(k_max,gamma,coeffs)
        %
        % Inputs
        %   n       Number of wavenumbers (integer)
        %   gamma   A parameter that determines frequency dependence
        %   a       If provided, then those coefficients are used to reproduce
        %           a GRF generated by a previous call to this function
        %
        % Outputs
        %   a               Fourier coefficients of GRF
        %   {V,Vp,Vpp,...}  Gaussian random field and it's 1st, 2nd, ..., nth
        %                   derivatives, respectively, in function form.
        
        % Generate new coefficients
        % or use existing coefficients
        if nargin < 3
            if nargin < 1, n = 50;     end
            if nargin < 2, gamma = 0;   end
            k = (0:n).';
            q = k/sqrt(n);
            sigma = sqrt(q.^gamma .* exp(-q.^2/2));
            sigma(1) = sigma(1) / sqrt(2);
            a = randn(length(q),2) .* [sigma sigma];
            a = a/norm(sigma([end:-1:2 1:end]))/sqrt(2);
        else
            k = (0:size(a,1)-1).';       % Wavenumber
            q = k/sqrt(numel(k)-1);      % Reduced wavenumber
        end
        
        % Assemble GRF and its derivatives
        varargout = cell(1,nargout-1);
        for nd = 0:nargout-2
%             if size(a,1) ~= 2, a = a.'; end
            aq = bsxfun(@times,a,q.^nd);
            switch mod(nd,4)
                case 0
                    f = @(x) sum(bsxfun(@times,  cos(q*x), aq(:,1) ),1) + ...
                             sum(bsxfun(@times,  sin(q*x), aq(:,2) ),1);
                case 1
                    f = @(x) sum(bsxfun(@times, -sin(q*x), aq(:,1) ),1) + ...
                             sum(bsxfun(@times,  cos(q*x), aq(:,2) ),1);
                case 2
                    f = @(x) sum(bsxfun(@times, -cos(q*x), aq(:,1) ),1) + ...
                             sum(bsxfun(@times, -sin(q*x), aq(:,2) ),1);
                case 3
                    f = @(x) sum(bsxfun(@times,  sin(q*x), aq(:,1) ),1) + ...
                             sum(bsxfun(@times, -cos(q*x), aq(:,2) ),1);
            end
            varargout{nd+1} = f;
        end
        
    end
    
    function [out] = build_template(c)
        
        varName = c(:,1).';
        val     = c{:,2};
        units   = c(:,3).';
        descrip = c(:,4).';
        
        out = table(c{:,2},'VariableNames',c(:,1).');
        out.Properties.VariableUnits = units;
        out.Properties.VariableDescriptions = descrip;
        
    end

end

end

function out = chop(a,b,fun)
    if nargin < 2, b = 0;     end
    if nargin < 3, fun = @le; end
    switch func2str(fun)
        case 'le'
            out = (a - b) <  1e-14;
        case 'ge'
            out = (a - b) > -1e-14;
    end
end

function out = structcat(varargin)
    out = varargin{1};
    for istruct = 2:nargin
        for fn = fieldnames(varargin{istruct}).'
            out.(fn{1}) = varargin{istruct}.(fn{1});
        end
    end
end

function varargout = build_potential(f,mv,mh)
    
    nd = min(length(f),nargout)-1;
    
    % prefactor = arrayfun(@(n) p.mv.^4 * p.mh^(-n),0:3,'Un',0);
    % potential = cellfun(@(grf,A) @(phi) A*grf(phi/mh),f,prefactor,'Un',0);
    % [V,Vp,Vpp] = deal(potential{1:3});
        
    for d = 0:nd
        grf = f{1+d};
        varargout{d+1} = @(phi) mv.^4 * mh^(-d) * reshape(grf(phi(:).'),size(phi));
    end
    
    for d = nd+1:max(nd,nargout-1)
        varargout{d+1} = function_handle.empty();
    end
    
end

%% Deprecated

% %     function [code] = gaussian_random_field_1D_py(a)
% %         
% %         k = 0:size(a,2)-1;            % Wavenumber
% %         q = k(:)'/sqrt(numel(k)-1);   % Reduced wavenumber
% %         
% %         code = sprintf('def V(phi):\n\treturn 0');
% %         return
% %         for kk = 1:length(k)
% %             code = strcat(code, ...
% %                 sprintf(' + cos(%.8f*phi) * %.8f',q(kk),a(1,kk)), ...
% %                 sprintf(' + sin(%.8f*phi) * %.8f',q(kk),a(2,kk)) );
% %         end
% %         
% %     end
    
% %         flag_testcase = false;
% %         if flag_testcase
% %             fname = 'C:\Users\Ross\OneDrive\Documents\Research\Eclipse Workspace\inf code\src\sample_coeffs_seed_m4.txt';
% %             data_raw = dlmread(fname,' ');
% %             kmax = ceil(length(data_raw)/2);
% %             a = zeros(2,kmax);
% %             a(1,:) = data_raw(kmax:-1:1);
% %             a(2,:) = data_raw(kmax:end);
% %             a(:,1) = a(:,1)/sqrt(2);
% %             p.Nafter = 1;
% %         end
% %         
% %         % Zoom in on the nearest minimum to approximate a quadratic potential
% %         flag_quadratic = false;
% %         if flag_quadratic
% %             V    = @(phi) (1/2)*p.mv^4*(phi/p.mh).^2;
% %             Vp   = @(phi) p.mv^4*(phi/p.mh)/p.mh;
% %             Vpp  = @(phi) p.mv^4/p.mh^2*ones(size(phi));
% %             Vppp = @(phi) zeros(size(phi));
% %             phistart = 10;
% %             p.Nafter = 1;
% %         end

% %         %% Stochastic Eternal Inflation with Adiabatic Regularization
% %         
% %         if false
% %             % what is m?
% %             % what to use for H*t?
% %             H = @(x) sqrt(V(x)/3); v = @(k,H,lna) k./H.*exp(-lna);
% %             m_H = @(x) m/H(x); n = @(x) sqrt(9/4 - m_H(x).^2);
% %             adiabaticRegSEIC = @(x) feval(@(v) sqrt(H.^2.*v.^3/(32*pi^2) .* (4*pi * abs(besselh(n,1,v))^2 - ...
% %                 (m_H(x)^2 + v^2)^(-7/2).*(8*m_H(x)^6 + 3*m_H(x)^4*(3 + 8*v^2) + ...
% %                 2*m_H(x)^2*v^2*(11 + 12*v^2) + 8*(v^4 + v^6)))),v(k,H,lna)) - 0.61*Delta_phi;
% %             
% %             % Find where stochastic eternal inflation starts/ends
% %             % Should always be eternal at phimax
% %             intv_ind  = find(diff(adiabaticRegSEIC(phi_N)));
% %             phibreak = arrayfun(@(i) fzero(adiabaticRegSEIC,phi_N(i:i+1)),intv_ind);
% %             
% %             if isempty(phibreak) || adiabaticRegSEIC(phiend)
% %                 obj.observables.NsinceEternal_AR = NaN; % Model is eternal up until phiend
% %             else
% %                 [~,imin] = min(abs(phibreak - phiend));
% %                 obj.observables.NsinceEternal_AR = integral(dlna_dphi,phibreak(imin),phistart);
% %             end
% %         end

%             phiedge = sqrt( 1 + 2*V0/Vpp0 + sqrt(-(4*V0 + Vpp0))*sqrt(-Vpp0)/Vpp0 );

% %     template = EternalInflationSimulator.build_template({...
% %             'timestamp',     nan('double'),     's',        'Timestamp of this experiment';...
% %             'status',        int8(-1),          '',         'Fate of this inflation scenario';...
% %             'Ntotal',       nan('single'),     '',         'Total number of e-foldings of inflation';...
% %             'Q',             nan('single'),     'M_pl',     'Amplitude of the scalar power spectrum at exit scale';...
% %             'r',             nan('single'),     '',         'Tensor-to-scalar ratio';...
% %             'n_s',           nan('single'),     '',         'Scalar spectral index';...
% %             'alpha',         nan('single'),     '',         'Running of the scalar spectral index';...
% %             'n_t',           nan('single'),     '',         'Tensor spectral index';...
% %             'dlgrho',        nan('single'),     '',         'Change in log(\rho) from exit scale to end of inflation';...
% %             'lgOk',          nan('single'),     '',         'log |\Omega - 1| (measure of flatness)';...
% %             'rho_Lambda',    nan('single'),     'M_pl^4',   'Vacuum energy of reheating minimum';...
% %             ...
% %             'NStochastic',   nan('single'),     '',         'Number of e-foldings of inflation during which stochastic eternal inflation criterion is valid';...
% %             'NSinceStoch',   nan('single'),     '',         'Number of e-foldings of inflation elapsed before starting point since stochastic eternal inflation criterion was valid';...
% %             'numStochEpochs',int8(-1),          '',         'Number of patches on the potential with stochastic eternal inflation';...
% %             ...
% %             'numTopolEpochs',int8(-1),          '',         'Number of patches on the potential with toplogical eternal inflation' });
    
% %         batch = repmat(struct(),1e3,1);
% %         
% %         for i_iter = 1:obj.parameters.n_iter
% %             
% %             obs = slowroll_with_tunneling(obj,0);
% %             
% %             % Time consuming
% %             for field = reshape(setdiff(fieldnames(obs),fieldnames(batch)),1,[])
% %                 [batch.(field{:})] = deal([]);
% %             end
% %             for field = reshape(setdiff(fieldnames(batch),fieldnames(obs)),1,[])
% %                 obs.(field{:}) = [];
% %             end
% %             if ~issorted(fieldnames(batch)), batch = orderfields(batch); end
% %             batch(mod(i_iter-1,1e3)+1) = orderfields(obs);
% %             
% %             if mod(i_iter,1e3) == 0
% %                 fprintf('Progress: %2.1f%%\n',(i_iter*100.0)/obj.parameters.n_iter);
% % %                 for field = reshape(setdiff(fieldnames(batch),fieldnames(obj.observables)),1,[])
% % %                     [obj.observables.(field{:})] = deal([]);
% % %                 end
% % %                 obs.observables(i_iter-1e3+(1:1e3)) = batch;
% %                 obj.observables{i_iter/1e3} = batch;
% %             end
% %             
% %             if mod(i_iter,1e5) == 0
% %                 save('ndim_slowroll_last_run','obj');
% %             end
% %             
% %         end

% %             %% Simulate slow roll
% %             
% %             % Generate potential as a Gaussian random field
% %             [a,f{1:4}] = obj.gaussian_random_field_1D(50,p.gamma);
% %             V = @(phi) p.mv^4*feval(@(grf) grf(phi/p.mh),f{1});
% %             
% %             [obs] = obj.slowroll(f);
% %             
% %             if n_tunnel == n_tunnel_max, break, end
             

% %     function [phitunnel] = check_false_vacuum_eternal_old(obj,a,potential,phi,mv)
% %         % Handle false vacuum tunneling.
% %         % Wrapper for Python-coded false-vacuum instanton calculation.
% %         
% %         p = obj.parameters;
% %         
% %         phistop   = phi(end);
% %         phitunnel = NaN;
% %         
% %         V = potential{1};
% %         
% % %         status = obj.compute_status(V,Vp,Vpp,phistop);
% % %         if status ~= 4, return; end
% % %         disp('Status = 4');
% %         
% %         % Find nearest neighbor minima to phistop
% %         phispace = phistop + linspace(-10*p.mh,10*p.mh,1001);
% %         [~,approx_mins] = findpeaks(-V(phispace));
% %         [~,imin] = min(abs(phispace(approx_mins)-phistop));
% %         nearest_mins = phispace(nonzeros(pos(approx_mins(imin+[-1 1]))));
% %         
% %         a = a.';
% %         
% %         p.V    = V;
% %         p.M_Pl = sqrt(8*pi);
% %         
% %         % 
% %         tunneling_rate = zeros(size(nearest_mins));
% %         for L_R = 1:numel(nearest_mins)
% %             
% %             phi_metaMin = phistop;
% %             phi_absMin  = nearest_mins(L_R);
% %             
% %             if V(phi_absMin) > V(phi_metaMin) || V(phi_absMin) < 0
% %                 continue
% %             end
% %             
% % %             a = py.list(a(:).');
% % %             try
% % %                 % Calculate Euclidean action for instanton solution
% % %                 S_E = py.InstantonWrapper.callInstanton(...
% % %                     a,mv,p.mh,phi_absMin,phi_metaMin);
% % %             catch err
% % %                 disp(err.message);
% % %                 continue
% % %             end
% %             
% %             %%
% %             p.phi_absMin    = +4.5e-1;
% %             p.phi_metaMin   = -4.5e-1;
% %             
% %             fvi = FalseVacuumInstanton(p);
% %             
% %             [R,Y,Rerr] = fvi.findProfile();
% %             S = fvi.findAction(R,Y);
% %             
% %             %%
% %             
% %             disp('Success');
% %             hbar = 1;
% %             preFactor = 1;
% %             
% %             % Calculate tunneling rate
% %             tunneling_rate(L_R) = preFactor*exp(-S_E/hbar);
% %             
% %         end
% %         
% %         % Choose "true" vacuum with the larger tunneling rate
% %         [tunneling_rate,imax] = max(tunneling_rate);
% %         
% %         datastruct.tunneling_rate = tunneling_rate;
% %         
% %         rateThreshold = 1e-5;
% %         if tunneling_rate < rateThreshold, return, end
% %         
% %         % Find the new value of phi post-tunneling
% %         phitunnel = findNewPhistart(nearest_mins(imax));
% %         
% %         datastruct.phitunnel = phitunnel;
% %         
% %     end
    