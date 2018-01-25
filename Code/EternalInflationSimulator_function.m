function EternalInflationSimulator_function(params,varargin)
% Constructor
% Inputs (varargin)
% 1) params     A structure containing some or all of the fields in
%               EternalInflationSimulator.parameters
% 2) {field1,val1,field2,val2,...}  Field and value pairs

%% Defaults

Mpl = sqrt(8*pi); % Planck mass

parameters = struct(...
    'n_iter',             1e2,...   % Number of iterations
    'mv',                 1e0,...   % Initial mass scale of the potential (multiple of Mpl)
    'mh',                 1,...     % Mass scale of the inflaton field (multiple of Mpl)
    'kmax',               30,...    % Largest wavenumber for GRF
    'gamma',              0,...     % Frequency dependence of GRF
    'Nafter',             55,...    % Number of e-folds between phiexit and phiend
    'lambdascreenmode',   true,...  % Throw out cases where rho_Lambda < 0?
    'fixLambda',          false,... % Condition on rho_Lambda ~= 0?
    'fixQ',               false,... % Condition on Q ~= 10^{-5}?
    'measure',            'B',...   % Measure on initial conditions
    'n_tunnel_max',       1,...     % Max number of tunneling events to simulate
    'outfile',            '');

if nargin > 0
    parameters = set_parameters(params,varargin{:});
end

main(parameters,Mpl)

end

function parameters = set_parameters(varargin)
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
                parameters.(fn{:}) = val.(fn{:});
            else
                error('parameters.n_iter must be a real, finite, positive integer.');
            end
            
        case 'mv'
            
            if isnumeric(val.mh) && isscalar(val.mh) && ...
                    isfinite(val.mh) && val.mh > 0 && isreal(val.mv)
                parameters.(fn{:}) = val.(fn{:});
            else
                error('parameters.mv must be a real, finite, positive number.');
            end
            
        case 'mh'
            
            if isnumeric(val.mh) && isscalar(val.mh) && ...
                    isfinite(val.mh) && val.mh > 0 && isreal(val.mh)
                parameters.(fn{:}) = val.(fn{:});
            else
                error('parameters.mh must be a real, finite, positive number.');
            end
            
        case 'kmax'
            if isnumeric(val.kmax) && isscalar(val.kmax) && ...
                    isfinite(val.kmax) && isreal(val.kmax)
                parameters.(fn{:}) = val.(fn{:});
            else
                error('parameters.kmax must be a real, finite number.');
            end
            
        case 'gamma'
            if isnumeric(val.gamma) && isscalar(val.gamma) && ...
                    isfinite(val.gamma) && isreal(val.gamma)
                parameters.(fn{:}) = val.(fn{:});
            else
                error('parameters.gamma must be a real, finite number.');
            end
            
        case 'Nafter'
            if isnumeric(val.Nafter) && isscalar(val.Nafter) && ...
                    isfinite(val.Nafter) && val.Nafter > 0 && ...
                    isreal(val.Nafter)
                parameters.(fn{:}) = val.(fn{:});
            else
                error('parameters.Nafter must be a finite, positive number.');
            end
            
        case 'lambdascreenmode'
            if islogical(val.lambdascreenmode) && isscalar(val.lambdascreenmode)
                parameters.(fn{:}) = val.(fn{:});
            else
                error('parameters.lambdascreenmode must be a logical.');
            end
            
        case 'fixLambda'
            if islogical(val.fixLambda) && isscalar(val.fixLambda)
                parameters.(fn{:}) = val.(fn{:});
            else
                error('parameters.fixLambda must be a logical.');
            end
            
        case 'fixQ'
            if islogical(val.fixQ) && isscalar(val.fixQ)
                parameters.(fn{:}) = val.(fn{:});
            else
                error('parameters.fixQ must be a logical.');
            end
            
        case 'measure'
            if ismember(upper(val.measure),{'A','B'})
                parameters.(fn{:}) = val.(fn{:});
            else
                error('parameters.measure must be ''A'' or ''B''');
            end
            
        case 'n_tunnel_max'
            if isnumeric(val.n_tunnel_max) && isscalar(val.n_tunnel_max) && ...
                    isfinite(val.n_tunnel_max) && val.n_tunnel_max > 0 && ...
                    mod(val.n_tunnel_max,1) == 0 && isreal(val.n_tunnel_max)
                parameters.(fn{:}) = val.(fn{:});
            else
                error('parameters.n_tunnel_max must be a real, finite, positive integer.');
            end
            
        case 'outfile'
            if ischar(val.outfile)
                parameters.outfile = val.outfile;
            else
                error('Output file name must be a string.');
            end
    end
end

end

function main(parameters,Mpl)
% Collect data from inflation simulations pertaining to the
% manifestation of false-vacuum, stochastic, and topological
% eternal inflation, as well as CMB observables.

    results_map = containers.Map({...
        'mv',...
        'status',...
        'Ntotal',...      
        'phitunnel',...    
        'log_tunneling_rate',...
        'Q',...          
        'r',...           
        'n_s',...        
        'alpha',...       
        'n_t',...        
        'dlgrho',...      
        'lgOk',...       
        'rho_Lambda',...  
        'NStochastic',...   
        'NSinceStoch',... 
        'numStochEpochs',...
        'numTopolEpochs'...
        },{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17});
    
p = parameters;

fid = fopen(p.outfile,'a');
if fid == -1
    fid = fopen(p.outfile,'w');
    fclose(fid);
    fid = fopen(p.outfile,'a');
end

tic

n_recycle = 4;
phistart_range = 8*p.mh*Mpl*...
    (-floor(n_recycle/2):floor(n_recycle/2)-1+mod(n_recycle,2));

disp(['outfile: ' p.outfile]);

% Run n_iter iterations of the inflation simulation
for i_iter = 1:p.n_iter
    
    data_out = nan(1,17);
    
    if mod(i_iter,1e3) == 0
        disp(num2str([i_iter toc]));
    end
    
    for goto = NaN % on break, goto record data
        
        data_out(1) = p.mv;
        
        %% Draw random potential and set initial conditions
        
        stopping_point = 0;
        
        % Only sample new potentials every 100 iterations
        if p.measure == 'B'
            
            % Generate potential function and
            % compute potential values at starting points
            if mod(i_iter,n_recycle) == 1
                [ak,f{1:3}] = gaussian_random_field_1D(p.kmax,p.gamma);
                Vstart   = p.mv^4*Mpl^4*f{1}(phistart_range/p.mh/Mpl);
                Vpstart  = p.mv^4*Mpl^3*f{2}(phistart_range/p.mh/Mpl)/p.mh;
                Vppstart = p.mv^4*Mpl^2*f{3}(phistart_range/p.mh/Mpl)/p.mh^2;
            end
            
            % Check if slow roll is valid at starting point
            ii = mod(i_iter-1,n_recycle)+1;
            if Vstart(ii) < 0
                data_out(2) = 1; break
            elseif (Vpstart(ii)/Vstart(ii)).^2/(16*pi/Mpl^2) > 1
                data_out(2) = 2; break
            elseif abs(Vppstart(ii)/Vstart(ii))/(8*pi/Mpl^2) > 1
                data_out(2) = 3; break
            else
                data_out(2) = 0;
            end
            phistart = phistart_range(ii);
            
        else
            
            [~,f{1:4}] = gaussian_random_field_1D(p.kmax,p.gamma);
            [data_init,phistart] = set_initial_conditions(f);
            data_out(2)          = data_init.status;
            
        end
        
        %% Simulate slow roll inflation
        
        stopping_point = 1;
        
        [status,mv,Ntotal,V,Vp,Vpp,phi,Vstop] = simulate_slowroll(parameters,Mpl,f,phistart,Vstart(ii),Vpstart(ii));
        
        data_out(1) = mv;
        data_out(2) = status;
        data_out(3) = Ntotal;
        
        %% Check for false-vacuum eternal inflation
        
        if isempty(Vstop), Vstop = V(phi(end)); end
        if status == 4 && Vstop > 0
            for i_tunnel = 1:p.n_tunnel_max
                
                % Look for an instanton tunneling solution
                [log_tunnel_rate,phitunnel,flag_hawking_moss] = check_false_vacuum_eternal(...
                    parameters,Mpl,f,phi,p.n_tunnel_max+1-i_tunnel,Vstop);
                
                if isnan(phitunnel)
                    break % No tunneling; move on
                elseif flag_hawking_moss
                    startAtMaximum = true;
                else
                    startAtMaximum = false;
                end
                
                stopping_point = 2;
                
                data_out(4) = phitunnel;
                data_out(5) = log_tunnel_rate;
                
                % If tunneling occurs, simulate slow roll from the new
                % starting point on other side of potential barrier
                [data_sr,V,Vp,Vpp,phi,Vstop] = simulate_slowroll(f,phitunnel,V(phitunnel));
                if isempty(Vstop), Vstop = V(phi(end)); end
                
                % Don't remember anything from before tunneling
                if false && SRStatus.compute_status(V,Vp,Vpp,phitunnel,Mpl) == 0 % Slow roll at exit point
                    if fast_tunnel
                        % Add pre-tunneling e-folds to total
                        data_sr.Ntotal = data_sr.Ntotal + Ntotal;
                    else
                        data_sr.Ntotal = Inf;
                        if data_sr.Ntotal < parameters.Nafter
                            phi(3) = phi_fv;
                        end
                    end
                end
                
            end
        end
        
        % Keep only the last bout of inflation
        data_out(1) = mv;
        data_out(2) = status;
        data_out(3) = Ntotal;
        
        if ~(Ntotal >= parameters.Nafter)
            break % Not enough inflation; skip to record data
        end
        
        %% Check for stochastic and topological eternal inflation
        
        stopping_point = 3;
        
        [numStochEpochs,NSinceStoch,NStochastic] = ...
            check_stochastic_eternal(parameters,Mpl,V,Vp,Vpp,phi,Ntotal);
        data_out(results_map('NSinceStoch'))    = NSinceStoch;
        data_out(results_map('numStochEpochs')) = numStochEpochs;
        data_out(results_map('NStochastic'))    = NStochastic;
        
        [numTopolEpochs] = check_topological_eternal(parameters,Mpl,V,Vp,Vpp,phi);
        data_out(results_map('numTopolEpochs')) = numTopolEpochs;
        
        %% Compute observables
        
        [~,f{4}] = gaussian_random_field_1D(p.kmax,p.gamma,ak);
        Vppp = build_potential(f,p.mv,p.mh,Mpl,3);
        
        Nbefore = Ntotal - p.Nafter; % e-folds before crossing
        observables = compute_observables(V,Vp,Vpp,Vppp,phi,Nbefore,Mpl);
        for fn = fieldnames(observables).'
            data_out(results_map(fn{1})) = observables.(fn{1});
        end
        
    end % for goto
    
    %% Record results
    
    data_out(isnan(data_out)) = 0;
    
    switch stopping_point
        %                 case 0
        %                     fprintf(fid,'%d\r\n',data_out(1));
        case 1
            fprintf(fid,'%d,%6.2f,%6.4f\r\n',data_out);
        case 2
            fprintf(fid,'%d,%6.2f,%6.4f,%6.4f,%6.4f\r\n',data_out);
        case 3
            fprintf(fid,'%d,%6.2f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%d,%d,%d,%d,%6.2f\r\n',data_out);
    end
    
end % for i_iter

fclose(fid);

end

%% Slowroll

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

startAtMaximum = false;

datastruct = obj.results_template;
datastruct.mv = p.mv;

[V,Vp,Vpp] = build_potential(f,p.mv,p.mh,obj.Mpl);

%% Initialize field value according to Measure A/B

switch upper(p.measure)
    
    case 'A'
        
        % Start at maximum (first slowroll region below max)
        if V(0) < 0
            datastruct.status = SRStatus(5);
            return
        end
        
        [phistart,status] = obj.find_phistart_A(0,V,Vp,Vpp,p.mh*obj.Mpl,obj.Mpl);
        if status == 4 && phistart ~= 0
            phistart = phistart + 1i; % Start at maximum
            startAtMaximum = true;
        end
        
    case 'B'
        
        phistart = 0;
        Vstart = V(phistart);
        if Vstart < 0
            datastruct.status = 1;
        elseif (Vp(phistart)/Vstart).^2/(16*pi/obj.Mpl^2) > 1
            datastruct.status = 2;
        elseif abs(Vpp(phistart)/Vstart)/(8*pi/obj.Mpl^2) > 1
            datastruct.status = 3;
        else
            datastruct.status = 0;
        end
        %                 datastruct.status = SRStatus.compute_status(V,Vp,Vpp,phistart,obj.Mpl);
        
end

end

function [status,mv,Ntotal,V,Vp,Vpp,phi,Vstop] = simulate_slowroll(parameters,Mpl,f,phistart,Vstart,Vpstart)
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

mv = parameters.mv;
mh = parameters.mh;

status = nan;
Ntotal = nan;

phi = nan(1,5);
phi(2) = phistart;

kappa = 8*pi/Mpl^2;

V = build_potential(f,mv,mh,Mpl);
Vp = []; Vpp = []; Vstop = [];

%         sgn = sign(Vpstart);
%         if sgn > 0
%             x = phistart + linspace(-4*mh*obj.Mpl,2*mh*obj.Mpl,101);
%         else
%             x = phistart + linspace(-2*mh*obj.Mpl,4*mh*obj.Mpl,101);
%         end
%         plot(x,V(x),'o'); hold on;
%         plot(phistart,V(phistart),'x'); hold off; pause(1)
%         V_interp = griddedInterpolant(x,V(x));

%% Simulate slow roll

shift_cases = [0];
if parameters.fixLambda, shift_cases = [shift_cases 1]; end
if parameters.fixQ,      shift_cases = [shift_cases 2]; end

for shifted_potential = shift_cases
    
    %% Modify potential to match observables
    
    switch shifted_potential
        
        case 1
            
            % If present-day cosmological constant is above some
            % threshold, subtract it from the potential and retry
            V_min = V(phistop);
            if abs(V_min) > 1e-20
                A = mv^4; grf = f{1};
                V = @(phi) A*grf(phi/mh*sqrt(kappa/8/pi)) - V_min;
            end
            
        case 2
            
            % Rescale mv so that Q matches observation
            Q_target = 2e-5;
            Q = sqrt(V(phiexit)/(150*(Vp(phiexit)./V(phiexit)).^2/2))/pi;
            mv = parameters.mv*sqrt(Q_target/Q);
            
            [V,Vp,Vpp] = build_potential(f,sqrt(Q_target/Q)*mv,mh,Mpl);
            
    end
    
    %% Identify phiend, phiexit, and phistop
    
    % Throw out models that settle into a negative minimum
    if parameters.lambdascreenmode
        % Just determine if V(phistop) < 0
        %                 phistop = EternalInflationSimulator.find_phistop_lite(phistart,V_interp,mh*obj.Mpl,Vstart,Vpstart);
        phistop = lambdascreen(phistart,V,mh*Mpl,Vstart,Vpstart);
        if isnan(phistop)
            status = 5;
            return
        end
    end
    
    % Only construct these when they're needed
    Vp  = build_potential(f,mv,mh,Mpl,1);
    Vpp = build_potential(f,mv,mh,Mpl,2);
    
    % Find the end of slow roll inflation
    [phiend,status] = find_phiend(phistart,V,Vp,Vpp,mh*Mpl,Mpl,Vstart,Vpstart);
    if status == 4 % Slow roll valid until minimum; no phiend
        phi(5) = phiend;
        return
    else
        phi(4) = phiend;
    end
    
    % Compute total number of slow roll e-folds
    dlna_dphi = @(phi) (-kappa*V(phi)./Vp(phi));
    switch upper(parameters.measure)
        case 'A'
            if startAtMaximum
                Ntotal = Inf;
            else
                Ntotal = integral(dlna_dphi,phistart,phiend);
            end
        case 'B'
            points = linspace(phistart,phiend,max(10,abs(phistart-phiend)/mv/Mpl/10));
            Ntotal_trapz = trapz(points,dlna_dphi(points));
            if  Ntotal_trapz > 0.7*parameters.Nafter && ...
                    Ntotal_trapz < 1.3*parameters.Nafter
                Ntotal = integral(dlna_dphi,phistart,phiend);
            else
                Ntotal = Ntotal_trapz;
            end
    end
    
    if Ntotal >= parameters.Nafter
        
        % Find the value of the field 55 e-folds
        % before the end of slow roll inflation
        dphi_dlna = @(N,phi) (-Vp(phi)./V(phi)/kappa).';
        Nbefore = Ntotal - parameters.Nafter;
        [~,phiexit] = ode45(dphi_dlna,[0 Nbefore Nbefore+1],phistart);
        phiexit = phiexit(2);
        phi(3) = phiexit;
        
    end
    
    % Find the nearest minimum after inflation ends
    [phistop,Vstop] = find_phistop(phiend,V,Vp,Vpp,mh*Mpl,false,true);
    phi(5) = phistop;
    
end

phimax = find_phistop(phistart,V,Vp,Vpp,mh*Mpl);
phi(1) = phimax;

if phistart == phiend
    return
end

% Move phistart to just shy of the local maximum
if strcmpi(parameters.measure,'A')
    delta = 1e-8*(mh*Mpl);
    if Vp(phistart) > 0, delta = -delta; end
    phistart = fminsearch(@(x) -V(x),phistart) + delta;
end

end

%% Checks for eternal inflation

function [log_tunnel_rate,phitunnel,flag_hawking_moss] = check_false_vacuum_eternal(parameters,Mpl,f,phi,n_tunnel_remaining,Vstop)
% Handle false vacuum tunneling.

if nargin < 5, n_tunnel_remaining = 1; end

mv = parameters.mv;
mh = parameters.mh;

phiscale = mh*Mpl;

%% Locate nearest minima

[V,Vp,Vpp] = build_potential(f,mv,mh,Mpl);
phistop = phi(end); % Location of false vacuum

if nargin < 5, Vstop = V(phistop); end

% Find minima around phistop
Vpstop = Vp(phistop);
near_minima(1) = find_next_minimum(phistop,V,Vp,Vpp,phiscale,-1,Vstop,Vpstop);
near_minima(2) = find_next_minimum(phistop,V,Vp,Vpp,phiscale,+1,Vstop,Vpstop);

%% Check for tunneling to nearby minima

new_phistart = nan(size(near_minima));
flag_hawking_moss = false(size(near_minima));
B = inf(size(near_minima));

for lr = find(~isnan(near_minima)) % Left and right neighbor basins
    
    %% Find barrier edge location
    
    phi_tol = abs(phistop-near_minima(lr))*1e-10;
    
    phimin = phistop;
    phimax = near_minima(lr);
    if phimin > phimax
        phitemp = phimin;
        phimin = phimax;
        phimax = phitemp;
    end
    phibar = 0.5*(phimin + phimax);
    phisep = abs(phimax-phimin);
    
    nbits = 5;
    ind = 1;
    while abs(phisep)*2^(1-ind) > abs(phi_tol)
        if mod(ind,nbits-1) == 1
            phisep = (phimax-phimin);
            phi_range = (phimin + phisep*(2^-nbits)):(phisep*(2^-nbits)):(phimax - phisep*(2^-nbits));
            fun_vals = Vstop-V(phi_range);
            ind = 1;
            ii = 2^(nbits-1);
        end
        if fun_vals(ii) > 0
            phimax = phi_range(ii);
            ii = ii - 2^(nbits-1-ind);
        else
            phimin = phi_range(ii);
            ii = ii + 2^(nbits-1-ind);
        end
        phibar = 0.5*(phimin+phimax);
        ind = ind + 1;
    end
    
    if (Vp(phibar)/Vstop).^2/(16*pi/Mpl^2) > 1
        status_bar = SRStatus.eps;
    elseif abs(Vpp(phibar)/Vstop)/(8*pi/Mpl^2) > 1
        status_bar = SRStatus.eta;
    else
        status_bar = SRStatus.ok;
    end
    
    %% Find the start of inflation in the neighboring basin
    
    if status_bar == 0 % No inflation at barrier edge
        phistart = phibar;
    else
        % Look for inflation further down the slope
        [phistart,status] = find_phistart_downhill(...
            phibar,V,Vp,Vpp,phiscale,Mpl);
        if status ~= 0, continue, end
    end
    
    %% Check if there is enough inflation
    
    % Find the end of slow roll inflation in neighbor basin
    [phiend,status] = find_phiend(phistart,V,Vp,Vpp,...
        phiscale,Mpl,[],[],false,n_tunnel_remaining > 1);
    
    % End of inflation does not occur before reacing the next
    % local minimum - cannot produce an observable universe in
    % that basin of the potential
    if n_tunnel_remaining == 1
        if status == 4, continue, end
    else
        if status ~= 4, continue, end
    end
    
    % Compute total number of slow roll e-folds
    dlna_dphi = @(phi) (-8*pi/Mpl^2*V(phi)./Vp(phi));
    points = linspace(phistart,phiend,max(10,abs(phistart-phiend)/mv/Mpl/10));
    Ntotal_trapz = trapz(points,dlna_dphi(points));
    if  Ntotal_trapz > 0.7*parameters.Nafter && ...
            Ntotal_trapz > 0*parameters.Nafter
        Ntotal = integral(dlna_dphi,phistart,phiend);
    else
        Ntotal = Ntotal_trapz;
    end
    
    % If we can't get enough e-foldings on the other side
    % the potential barrier, don't bother computing the
    % tunneling rate
    if Ntotal < parameters.Nafter
        continue
    end
    
    continue
    
    %% Do full instanton calculation
    
    % Define potential with rescaled mv = 1
    [V_rescale,Vp_rescale,Vpp_rescale] = build_potential(f,1,mh,Mpl);
    
    xtol         = 1e-4;
    phitol       = 1e-4;
    thinCutoff   = 1e-2;
    
    % Initialize instanton solver
    fvi = FalseVacuumInstanton(...
        'V',            V_rescale,...
        'dV',           Vp_rescale,...
        'd2V',          Vpp_rescale,...
        'M_Pl',         Mpl,...
        'phi_metaMin',  phistop,...
        'phi_absMin',   near_minima(lr));
    
    try % Solve for instanton profile
        [R,Y,~] = fvi.find_profile([],xtol,phitol,thinCutoff);
    catch me
        switch me.identifier
            case 'FalseVacuumInstanton:StableFalseVacuum'
                continue % No tunneling
            otherwise
                rethrow(me);
        end
    end
    
    if isscalar(R)
        flag_hawking_moss(lr) = true;
    end
    
    new_phistart(lr) = Y(1,1); % Field value at center of bubble
    
    % Get tunneling suppression rate B = -log(\lambda)
    % using appropriately scaled radial coordinate
    B(lr) = fvi.find_tunneling_suppression(R/mv^2,Y);
    
end

%% Choose bubble with the larger tunneling rate

[log_tunnel_rate,imax] = max(-B);

if all(isnan(new_phistart))
    phitunnel = NaN;
    return
end

%% Find the new value of phi after tunneling

kappa = 8*pi/Mpl^2;
if log_tunnel_rate >= log(9/4/pi) + 2*log(kappa/3*V_rescale(phistop)*mv^4);
    phitunnel = new_phistart(imax);
    flag_hawking_moss = flag_hawking_moss(imax);
else
    phitunnel = NaN;
    flag_hawking_moss = false;
end

end
    
function [numStochEpochs,NSinceStoch,NStochastic] = check_stochastic_eternal(parameters,Mpl,V,Vp,Vpp,phi,Ntotal)

numStochEpochs = nan;
NSinceStoch    = nan;
NStochastic    = nan;

p = parameters;

if isnumeric(phi), phi = num2cell(phi); end

[phimax,phistart,phiexit,phiend] = deal(phi{1:4});

if isnan(phimax)
    phimax = find_phimax(phistart,V,Vp,Vpp,p.mh*Mpl);
end

dlna_dphi = @(phi) (-V(phi)./Vp(phi));
dphi_dlna = @(N,phi) (-Vp(phi)./V(phi)).';

kappa = 8*pi/Mpl^2;
% \kappa = 1
slowroll = @(x) chop((Vp(x)./V(x)).^2/2/kappa,1) & chop(abs(Vpp(x)./V(x)/kappa),1);

%% Stochastic Eternal inflation

% Here \kappa = 1
stochasticEIC = @(x) (kappa*V(x)).^(3/2) > 2*pi*sqrt(3)*(0.607)*abs(Vp(x)); % stochasticEIC > 0 -> Eternal

% Find the value of the field 55 e-folds
% before the end of slow roll inflation

lna_steps = 0:floor(Ntotal);
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

%         disp(num2str(length(phibreak_sei)));

%         disp([num2str(V(phistart)^(3/2)) ' ' num2str(2*pi*sqrt(3)*(0.607)*abs(Vp(phistart)))]);

if isempty(phibreak_sei)
    phibreak_sei = phimax;
    off2on_sei = 0;
end

% Second derivative check for stochastic eternal inflation near the
% maximum. Assume potential is locally mirror symmetric about maximum,
% or we condition on the field falling only toward phistart.
phidev = phimax-phibreak_sei(1); Hmax = sqrt(V(phimax)/3);
if abs(erf(abs(phidev/(Hmax/2/pi))/sqrt(2))) > exp(-3)
    numStochEpochs = uint8(1 + nnz(off2on_sei));
else
    stochasticEIC = @(x) stochasticEIC(x) & ...
        ~(min(phimax,phibreak_sei(1)) <= x && x <= max(phimax,phibreak_sei(1)));
    phibreak_sei(1) = []; off2on_sei(1) = [];
    numStochEpochs = uint8(nnz(off2on_sei));
end

% Compute # of e-folds after past SEI breakdown and before phiexit
if stochasticEIC(phiexit)
    NSinceStoch = 0; % Eternal at phiexit
elseif any((phibreak_sei(~off2on_sei)-phiexit)*Vp(phiexit) > 0)
    phib = phibreak_sei((phibreak_sei-phiexit)*Vp(phiexit) > 0 & ~off2on_sei);
    [~,imin] = min(abs(phib-phiexit));
    NSinceStoch = integral(@(x) dlna_dphi(x).*slowroll(x),...
        phib(imin),phiexit);
else
    NSinceStoch = NaN; % Never eternal before phiexit
end

% Compute # of e-folds eternal after phistart and before phiend
if strcmpi(p.measure,'a') && stochasticEIC(phistart)
    NStochastic = Inf;
else
    NStochastic = nan;
    %             datastruct.NStochastic = integral(@(x) ...
    %                 dlna_dphi(x).*stochasticEIC(x).*slowroll(x),phistart,phiend);
end

end
    
function [numTopolEpochs] = check_topological_eternal(parameters,Mpl,V,Vp,Vpp,phi)

numTopolEpochs = nan;

phimax   = phi(1);
phistart = phi(2);

%% Topological Eternal Inflation

kappa = 8*pi/Mpl^2;

% Determine whether quantum fluctuations could result in at least
% one Hubble volume descending toward a different minimum of V(phi)
if strcmpi(parameters.measure,'a')
    checkTopological = true;
else
    Vstart = V(phistart);
    if Vstart > 0
        H = sqrt(kappa*Vstart/3);
        dphi = H/2/pi;
        Dphi = -Vp(phistart)/(3*H^2);
        try
            checkTopological = 1/2*erfc(abs(phistart+Dphi-phimax)/(sqrt(2)*dphi)) > exp(-3);
        catch me
            disp('');
        end
    else
        checkTopological = false;
    end
end

if checkTopological
    
    phimax = phi{1};
    
    % Find value of phi at the domain wall boundary
    phiedge_eps = fzero(@(phi) Vp(phi)./V(phi)/2/kappa - 1,phimax);
    phiedge_eta = fzero(@(phi) abs(Vpp(phi)./V(phi)/kappa) - 1,phimax);
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
        numTopolEpochs = true;
    end
    
end

end

%% Subroutines

function [phiend,status] = find_phiend(phistart,V,Vp,Vpp,phiscale,Mpl,Vstart,Vpstart,lambdascreenmode,precisephistop)
% Find the value of phi at the end of slow roll inflation

if nargin < 7 || isempty(Vstart), Vstart = V(phistart);     end
if nargin < 8 || isempty(Vpstart), Vpstart = Vp(phistart);   end
if nargin < 9, lambdascreenmode = false; end
if nargin < 10, precisephistop = true; end

status = 0;

sgn    = sign(Vpstart);
phimin = phiscale^2 * abs(Vpstart./Vstart);
dphi   = -0.001*sgn*max(0.01*phiscale,min(phimin,phiscale));            %%% Tunable

% Take steps until passed breakdown of slow roll
step = 10^(1/16);                                                       %%% Tunable
ind = 1;
batch = 40;
phi = phistart*ones(1,batch);
while status == 0
    if mod(ind,batch) == 1
        phi_last = phi(end);
        phi = phi_last + cumsum(dphi*step.^(ind-1:ind+batch-2));
        Vend   = V(phi);
        Vpend  = Vp(phi);
        Vppend = Vpp(phi);
    end
    ii = mod(ind-1,batch)+1;
    if Vend(ii) < 0
        status = 1;
    elseif (Vpend(ii)/Vend(ii)).^2/(16*pi/Mpl^2) > 1
        status = 2;
    elseif abs(Vppend(ii)/Vend(ii))/(8*pi/Mpl^2) > 1
        status = 3;
    end
    if sgn*Vpend(ii) < 0
        status = 4;
    end
    ind = ind + 1;
end
phiend = phi;

if lambdascreenmode
    if status == 1
        phiend = phi(ii);
        return
    end
    % Continue search until minimum reached or negative potential
    while true
        if mod(ind,batch) == 1
            phi_last = phi(end);
            phi = phi_last + cumsum(dphi*step.^(ind-1:ind+batch-2));
            Vend   = V(phi);
            Vpend  = Vp(phi);
            Vppend = Vpp(phi);
        end
        ii = mod(ind-1,batch)+1;
        %                 plot(x,V(x)); hold on;
        %                 plot(phi(ii),V(phi(ii)),'o'); hold off;
        %                 pause(0.1);
        if Vend(ii) < 0
            phiend = phi(ii);
            status = status + 5;
            return;
        elseif sgn*Vpend(ii) < 0
            break;
        end
        ind = ind + 1;
    end
end

% Choose function to zero based on status
switch status
    case 1, fun = V;                            % V = 0
    case 2, fun = @(x) (Vp(x)./V(x)).^2/2 - 1;  % eps = 1
    case 3, fun = @(x) abs(Vpp(x)./V(x)) - 1;   % eta = 1
    case 4, fun = Vp;                           % Vp = 0
end

phimax = phiend(ii);
if ii > 1
    phimin = phiend(ii-1);
else
    phimin = phi_last;
end
if status == 1 || (status == 4 && sgn == 1)
    phitmp = phimin;
    phimin = phimax;
    phimax = phitmp;
end

phisep = (phimax-phimin);
phiend = 0.5*(phimin+phimax);

if ~precisephistop && status == 4
    return
end

nbits = 5;
ind = 1;
while abs(phisep)*2^(1-ind) > abs(dphi)
    if mod(ind,nbits-1) == 1
        phisep = (phimax-phimin);
        phi = (phimin + phisep*(2^-nbits)):(phisep*(2^-nbits)):(phimax - phisep*(2^-nbits));
        fun_vals = fun(phi);
        ind = 1;
        ii = 2^(nbits-1);
    end
    if fun_vals(ii) > 0
        phimax = phi(ii);
        ii = ii - 2^(nbits-1-ind);
    else
        phimin = phi(ii);
        ii = ii + 2^(nbits-1-ind);
    end
    phiend = 0.5*(phimin+phimax);
    ind = ind + 1;
end

end

function [phistop,Vstop] = find_phistop(phiend,V,Vp,~,phiscale,lambdascreenmode,precisionmode)
% Find the value of phi at the next local minimum

if nargin < 6, lambdascreenmode = false; end
if nargin < 7, precisionmode = true;     end

sgn    = sign(Vp(phiend));
phimin = phiscale^2 * abs(Vp(phiend)./V(phiend));
dphi   = -0.001*sgn*max(0.01*phiscale,min(phimin,phiscale));            %%% Tunable

%         x = linspace(-5*phiscale,5*phiscale,1001);

Vstop = [];

% Take steps until passed local min
step = 10^(1/16);                                                       %%% Tunable
step = 1e2;
ind = 1;
batch = 30;
phistop = phiend*ones(1,batch);                                                      %%% Tunable
% cumsum(dphi*step.^(ind-1:ind+batch-2))
while true
    if mod(ind,batch) == 1
        phistop_last = phistop(end);
        phistop = phistop_last + cumsum(dphi*step*(ind-1:ind+batch-2));
        Vpstop  = Vp(phistop);
        if lambdascreenmode
            Vstop   = V(phistop);
        end
    end
    ii = mod(ind-1,batch)+1;
    %             plot(x,V(x)); hold on;
    %             plot(phistop(ii),V(phistop(ii)),'o'); hold off;
    %             pause(0.1);
    if lambdascreenmode && Vstop(ii) < 0
        phistop = nan;
        return
    end
    if sgn*Vpstop(ii) < 0
        break
    end
    ind = ind + 1;
end

phimin = phistop(ii);
if ii > 1
    phimax = phistop(ii-1);
else
    phimax = phistop_last;
end

if precisionmode
    
    phistop = 0.5*(phimin + phimax);
    while abs(phimax-phimin) > abs(dphi)
        Vpstop = Vp(phistop);
        if sgn*Vpstop > 0
            phimax = phistop;
        else
            phimin = phistop;
        end
        phistop = 0.5*(phimin + phimax);
    end
    
else
    
    phistop = [phimin,phimax];
    
end

end

function [phistop] = lambdascreen(phiend,V,phiscale,Vstart,Vpstart)
% Find the value of phi at the next local minimum

if nargin < 4, Vstart = V(phistart);     end
if nargin < 5, Vpstart = Vp(phistart);   end

sgn    = sign(Vpstart);
phimin = phiscale^2 * abs(Vpstart./Vstart);
dphi   = -0.001*sgn*max(0.01*phiscale,min(phimin,phiscale));            %%% Tunable

% Take steps until passed local min
step = 1e2;
ind = 1;
batch = 30;
phistop = phiend*ones(1,batch);
Vstop = inf(1,batch);
while true
    if mod(ind,batch) == 1
        phistop_last = phistop(end);
        Vstop_last = Vstop(end);
        phistop = phistop_last + cumsum(dphi*step*(ind-1:ind+batch-2));
        Vstop = V(phistop);
    end
    ii = mod(ind-1,batch)+1;
    if Vstop(ii) < 0
        phistop = nan;
        return
    end
    if ii > 1 && Vstop(ii) > Vstop(ii-1)
        break
    elseif ii == 1 && Vstop(1) > Vstop_last
        break
    end
    ind = ind + 1;
end

phimin = phistop(ii);
if ii > 1
    phimax = phistop(ii-1);
else
    phimax = phistop_last;
end
phistop = 0.5*(phimin+phimax);

end

function [phipeak] = find_phimax(phi0,V,Vp,~,phiscale)
% Find the value of phi at the previous local maximum

sgn    = sign(Vp(phi0));
phimin = phiscale^2 * abs(Vp(phi0)./V(phi0));
dphi   = -0.001*sgn*max(0.01*phiscale,min(phimin,phiscale));            %%% Tunable

% Take steps until passed local max
step = 1e2;
ind = 1;
batch = 30;
phi = phi0*ones(1,batch);                                                      %%% Tunable
while true
    if mod(ind,batch) == 1
        phi_last = phi(end);
        phi = phi_last - cumsum(dphi*step*(ind-1:ind+batch-2));
        Vpstop  = Vp(phi);
    end
    ii = mod(ind-1,batch)+1;
    if sgn*Vpstop(ii) < 0
        break
    end
    ind = ind + 1;
end

phimin = phi(ii);
if ii > 1
    phimax = phi(ii-1);
else
    phimax = phi_last;
end
phipeak = 0.5*(phimin + phimax);
while abs(phimax-phimin) > abs(dphi)
    if sgn*Vp(phipeak) > 0
        phimax = phipeak;
    else
        phimin = phipeak;
    end
    phipeak = 0.5*(phimin + phimax);
end

end

function [phistart,status] = find_phistart_downhill(phiinit,V,Vp,Vpp,phiscale,Mpl,status0)
% Find the value of phi at the end of slow roll inflation

sgn    = sign(Vp(phiinit));
phimin = phiscale^2 * abs(Vp(phiinit)./V(phiinit));
dphi   = -0.001*sgn*max(0.01*phiscale,min(phimin,phiscale));

status = 2;

step = 10^(1/16);
ind = 1;
batch = 40;
phi = phiinit*ones(1,batch);
while status ~= 0
    status_last = status;
    if mod(ind,batch) == 1
        phi_last = phi(end);
        phi = phi_last + cumsum(dphi*step.^(ind-1:ind+batch-2));
        Vend   = V(phi);
        Vpend  = Vp(phi);
        Vppend = Vpp(phi);
    end
    ii = mod(ind-1,batch)+1;
    if (Vpend(ii)/Vend(ii)).^2/(16*pi/Mpl^2) > 1
        status = 2;
    elseif abs(Vppend(ii)/Vend(ii))/(8*pi/Mpl^2) > 1
        status = 3;
    else
        status = 0;
    end
    if sgn*Vpend(ii) < 0
        status = 4;
        break
    end
    ind = ind + 1;
end

if status ~= 0
    phistart = nan;
    return
end

%% Find the precise starting point binary search

% Choose function to zero based on status
switch status_last
    case 1, fun = V;                            % V = 0
    case 2, fun = @(x) (Vp(x)./V(x)).^2/2 - 1;  % eps = 1
    case 3, fun = @(x) abs(Vpp(x)./V(x)) - 1;   % eta = 1
    case 4, fun = Vp;                           % Vp = 0
end

phimax = phi(ii);
if ii > 1
    phimin = phi(ii-1);
else
    phimin = phi_last;
end

if phimin > phimax
    phitemp = phimin;
    phimin = phimax;
    phimax = phitemp;
end

phistart = 0.5*(phimin + phimax);
while abs(phimax-phimin) > abs(dphi)
    if fun(phistart) > 0
        phimax = phistart;
    else
        phimin = phistart;
    end
    phistart = 0.5*(phimin + phimax);
end

end

function [phinextmin] = find_next_minimum(phistop,V,Vp,Vpp,phiscale,search_direction,Vstop,Vpstop)
% Find the value of phi at the next local minimum

if nargin < 6, search_direction = 1; end

phimin = phiscale^2 * abs(Vpstop./Vstop);
dphi   = 0.001*search_direction*max(0.01*phiscale,min(phimin,phiscale));

%         x = linspace(-5*phiscale,5*phiscale,1001);

% Take steps until passed local max
flag_crossed_maximum = false;
step = 3e3;
ind = 1;
batch = 25;
phinextmin = phistop*ones(1,batch);
%         dphi*step.^(ind-1:ind+batch-2)
while true
    if mod(ind,batch) == 1
        phinextmin_last = phinextmin(end);
        phinextmin = phinextmin_last + cumsum(dphi*step*(ind-1:ind+batch-2));
        if flag_crossed_maximum
            Vstop = V(phinextmin);
        end
        Vpstop  = Vp(phinextmin);
        Vppstop = Vpp(phinextmin);
    end
    ii = mod(ind-1,batch)+1;
    %             plot(x,V(x)); hold on;
    %             plot(phinextmin(ii),V(phinextmin(ii)),'o'); hold off;
    %             pause(0.1);
    if flag_crossed_maximum
        if Vstop(ii) < 0
            phinextmin = nan;
            return
        elseif search_direction*Vpstop(ii) > 0
            break
        end
    elseif search_direction*Vpstop(ii) < 0 && Vppstop(ii) < 0
        flag_crossed_maximum = true;
        if ii < batch
            Vstop = [zeros(1,ii) V(phinextmin(ii+1:end))];
        end
    end
    ind = ind + 1;
end

phimin = phinextmin(ii);
if ii > 1
    phimax = phinextmin(ii-1);
else
    phimax = phinextmin_last;
end

if phimin > phimax
    phitemp = phimin;
    phimin = phimax;
    phimax = phitemp;
end
phinextmin = 0.5*(phimin + phimax);
phisep = abs(phimax-phimin);

nbits = 5;
ind = 1;
while abs(phisep)*2^(1-ind) > abs(dphi)
    if mod(ind,nbits-1) == 1
        phisep = (phimax-phimin);
        phi = (phimin + phisep*(2^-nbits)):(phisep*(2^-nbits)):(phimax - phisep*(2^-nbits));
        fun_vals = Vp(phi);
        ind = 1;
        ii = 2^(nbits-1);
    end
    if fun_vals(ii) > 0
        phimax = phi(ii);
        ii = ii - 2^(nbits-1-ind);
    else
        phimin = phi(ii);
        ii = ii + 2^(nbits-1-ind);
    end
    phinextmin = 0.5*(phimin+phimax);
    ind = ind + 1;
end

end

function [phistart,status] = find_phistart_A(phi0,V,Vp,Vpp,phiscale,Mpl)
% Find the value of phi at the end of slow roll inflation

% Screen negative starting points?

orig_status = SRStatus.compute_status(V,Vp,Vpp,phi0,Mpl);
status = orig_status;

sgn    = sign(Vp(phi0));
phimin = phiscale^2 * abs(Vp(phi0)./V(phi0));
dphi   = -0.001*sgn*max(0.01*phiscale,min(phimin,phiscale));

sr_at_0 = status.isok(); % Is SRA valid at phistart?
if sr_at_0
    search_direction = -1; % If inflating, look up
else
    search_direction = 1;  % If not inflating, look down
end

% Take steps until passed breakdown of slow roll
step = 10^(1/16);
while (sr_at_0 && status.isok()) || (~sr_at_0 && status.isok())
    phistart = phi0 + search_direction*dphi;
    status = SRStatus.compute_status(V,Vp,Vpp,phistart,Mpl);
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
    case SRStatus.rho, fun = V;                                             % V = 0
    case SRStatus.eps, fun = @(x) 1 - (Vp(x)./V(x)).^2/(16*pi/obj.M_pl^2);  % eps = 1
    case SRStatus.eta, fun = @(x) 1 - abs(Vpp(x)./V(x))/(8*pi/obj.M_pl^2);  % eta = 1
    case SRStatus.localmin, fun = Vp;                                       % Vp = 0
end

% Find zero of function to obtain true phistart
phistart = fzero(fun,[phistart phi0]);

end

function [status] = compute_status(V,Vp,Vpp,phi,Mpl)
if V(phi) < 0
    status = 1;
elseif (Vp(phi)/V(phi)).^2/(16*pi/Mpl^2) > 1
    status = 2;
elseif abs(Vpp(phi)/V(phi))/(8*pi/Mpl^2) > 1
    status = 3;
else
    status = 0;
end
end

function [observables] = compute_observables(V,Vp,Vpp,Vppp,phi,Nbefore,Mpl)
% Inputs
%   potential   Cell array of function handles
%                   {V, V', V'', V'''}
%   phi         Monotonic array of field values
%                   [phistart, phiexit, phiend, phistop]
%   Nbefore     Number of e-foldings between phistart and phiexit
%
% Outputs
%   observables Structure of observables

%         phistart = phi(2);
phiexit  = phi(3);
%         phiend   = phi(4);
%         phistop  = phi(5);

V_exit = V(phiexit);

kappa = 8*pi/Mpl^2;

% Dimensionless slow roll parameters at horizon exit scale
eps = (Vp(phiexit)./V_exit).^2/(2*kappa);
eta = Vpp(phiexit)./V_exit/kappa;
xi2 = Vp(phiexit).*Vppp(phiexit)./V_exit^2/kappa^2;

observables.Q           = sqrt(V_exit/(150*eps))/pi;
observables.r           = 16*eps;
observables.n_s         = 1-6*eps+2*eta;
observables.alpha       = 16*eps*eta - 24*eps^2 - 2*xi2;
observables.n_t         = -2*eps;
observables.dlgrho      = log10(V_exit/V(phi(4)));
observables.lgOk        = log10(V(phi(3))/V_exit) - Nbefore*2/log(10);
observables.rho_Lambda  = V(phi(5));

end

%% Utilities

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
    aq = bsxfun(@times,a,q.^nd);
    switch mod(nd,4)
        case 0
            f = @(x) sum(bsxfun(@times,  cos(q*x), aq(:,1) ) + ...
                bsxfun(@times,  sin(q*x), aq(:,2) ),1);
        case 1
            f = @(x) sum(bsxfun(@times, -sin(q*x), aq(:,1) ) + ...
                bsxfun(@times,  cos(q*x), aq(:,2) ),1);
        case 2
            f = @(x) sum(bsxfun(@times, -cos(q*x), aq(:,1) ) + ...
                bsxfun(@times, -sin(q*x), aq(:,2) ),1);
        case 3
            f = @(x) sum(bsxfun(@times,  sin(q*x), aq(:,1) ) + ...
                bsxfun(@times, -cos(q*x), aq(:,2) ),1);
    end
    varargout{nd+1} = f;
end

end

function [a,varargout] = gaussian_random_field_1D_single(n,gamma,a)
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
    q = single(k/sqrt(n));
    sigma = sqrt(q.^gamma .* exp(-q.^2/2));
    sigma(1) = sigma(1) / sqrt(2);
    a = randn(length(q),2) .* [sigma sigma];
    a = single(a/norm(sigma([end:-1:2 1:end]))/sqrt(2));
else
    k = (0:size(a,1)-1).';       % Wavenumber
    q = single(k/sqrt(numel(k)-1));      % Reduced wavenumber
end

% Assemble GRF and its derivatives
varargout = cell(1,nargout-1);
for nd = 0:nargout-2
    aq = bsxfun(@times,a,q.^nd);
    switch mod(nd,4)
        case 0
            f = @(x) sum(bsxfun(@times,  cos(q*x), aq(:,1) ) + ...
                bsxfun(@times,  sin(q*x), aq(:,2) ),1);
        case 1
            f = @(x) sum(bsxfun(@times, -sin(q*x), aq(:,1) ) + ...
                bsxfun(@times,  cos(q*x), aq(:,2) ),1);
        case 2
            f = @(x) sum(bsxfun(@times, -cos(q*x), aq(:,1) ) + ...
                bsxfun(@times, -sin(q*x), aq(:,2) ),1);
        case 3
            f = @(x) sum(bsxfun(@times,  sin(q*x), aq(:,1) ) + ...
                bsxfun(@times, -cos(q*x), aq(:,2) ),1);
    end
    varargout{nd+1} = f;
end

end

function [a,varargout] = gaussian_random_field_1D_at_x(n,gamma,phi)
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
if nargin < 1, n = 50;     end
if nargin < 2, gamma = 0;   end
k = (0:n).';
q = k/sqrt(n);
sigma = sqrt(q.^gamma .* exp(-q.^2/2));
sigma(1) = sigma(1) / sqrt(2);
a = randn(length(q),2) .* [sigma sigma];
a = a/norm(sigma([end:-1:2 1:end]))/sqrt(2);

% Assemble GRF and its derivatives
varargout = cell(1,nargout-1);
qphi = q*phi;
for nd = 0:nargout-2
    %             if size(a,1) ~= 2, a = a.'; end
    if nd == 0
        aq = a;
    else
        aq = bsxfun(@times,aq,q);
    end
    switch mod(nd,4)
        case 0
            f = sum(bsxfun(@times,  cos(qphi), aq(:,1) ),1) + ...
                sum(bsxfun(@times,  sin(qphi), aq(:,2) ),1);
        case 1
            f = sum(bsxfun(@times, -sin(qphi), aq(:,1) ),1) + ...
                sum(bsxfun(@times,  cos(qphi), aq(:,2) ),1);
        case 2
            f = sum(bsxfun(@times, -cos(qphi), aq(:,1) ),1) + ...
                sum(bsxfun(@times, -sin(qphi), aq(:,2) ),1);
        case 3
            f = sum(bsxfun(@times,  sin(qphi), aq(:,1) ),1) + ...
                sum(bsxfun(@times, -cos(qphi), aq(:,2) ),1);
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

function varargout = build_potential(f,mv,mh,Mpl,d)
    
    nd = min(length(f),nargout)-1;
    
    mv = mv*Mpl;
    mh = mh*Mpl;
    
    % prefactor = arrayfun(@(n) p.mv.^4 * p.mh^(-n),0:3,'Un',0);
    % potential = cellfun(@(grf,A) @(phi) A*grf(phi/mh),f,prefactor,'Un',0);
    % [V,Vp,Vpp] = deal(potential{1:3});
    
    mv4 = mv.^4;
    
    if nargin < 5
        for d = 0:nd
            grf = f{1+d};
            mhmd = mh^(-d);
            varargout{d+1} = @(phi) mv4 * mhmd * grf(phi/mh); % reshape(grf(phi(:).'),size(phi));
        end
    else
        grf = f{1+d};
        mhmd = mh^(-d);
        varargout{1} = @(phi) mv4 * mhmd * grf(phi/mh); % reshape(grf(phi(:).'),size(phi));
    end
    
    for d = nd+1:max(nd,nargout-1)
        varargout{d+1} = function_handle.empty();
    end
    
end

function x = binarysearch(fun,xmin,xmax,dx,swap_condition)

if nargin < 5, swap_condition = false; end

if swap_condition
    xtmp = xmin;
    xmin = xmax;
    xmax = xtmp;
end

x = 0.5*(xmin + xmax);

while abs(xmax-xmin) > abs(dx)
    if fun(x) > 0
        xmax = x;
    else
        xmin = x;
    end
    x = 0.5*(xmin + xmax);
end

end
