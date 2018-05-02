classdef EternalInflationSimulator < handle
% Running EternalInflationSimulator.main() performs a set of inflation
% simulations and records data relevant to the presence of eteranl
% inflation. At each iteration, the code generates a potential as a
% Gaussian random field and initializes the inflaton field according to
% some specified measure.

%% Properties, constructor, and set/get methods

properties
    
    parameters = struct(...
        'seed',               NaN,...   % Seed for random number generator (integer)
        'n_iter',             1e2,...   % # of iterations
        'mv',                 1e0,...   % Initial mass scale of the potential (multiple of m_Pl)
        'mh',                 1,...     % Mass scale of the inflaton field (multiple of m_Pl)
        'kmax',               30,...    % Largest wavenumber for GRF
        'gamma',              0,...     % Frequency dependence of GRF
        'Nafter',             55,...    % Number of e-folds between phiexit and phiend
        'lambdascreenmode',   true,...  % Throw out cases where rho_Lambda < 0?
        'fixQ',               false,... % Condition on Q ~= 2*10^{-5}?
        'measure',            'B',...   % Measure on initial conditions
        'n_tunnel_max',       1,...     % Max number of tunneling events to simulate
        'n_recycle',          4,...     % # of times to reuse same V(phi) with different phi0 value
        'rho_Lambda_thres',   0,...     % Threshold below which vacuum energy is considered "small"
        'outfile',            ''...     % Path to output text file
        );
    
    plotFlag = 0;
    plotHandle = [];
    
end

properties (Constant)
    
    % Always hbar = c = 1
    m_Pl = sqrt(8*pi); % Planck mass
    M_Pl = 1;  % Reduced Planck mass
    
end

%% Main method

methods (Access = public)
    
    function main(obj)
        % Collect data from inflation simulations pertaining to the
        % manifestation of false-vacuum, stochastic, and topological
        % eternal inflation, as well as CMB observables.
        
        p = obj.parameters;
        
        if p.n_tunnel_max < 0
            p.flag_screen_basin = true;
            p.n_tunnel_max = -1*p.n_tunnel_max;
        else
            p.flag_screen_basin = false;
        end
        
        % Create output text file, write metadata header
        fid = fopen(p.outfile,'wt');
        fprintf(fid,'%E,%.4G,%.4G,%.4G,%d,%.2f,%s,%d,%d,%.4G,%d,%.3f,%d,%d\r\n',...
            p.n_iter,p.mv,p.mh,obj.m_Pl,p.kmax,p.gamma,p.measure,p.n_tunnel_max,...
            p.lambdascreenmode,p.rho_Lambda_thres,p.fixQ,p.Nafter,p.seed,p.n_recycle);
        fclose(fid);
        
        % Seed the random number generator
        if ~isnan(p.seed), rng(p.seed); end
        
        % Mass scales in natural units
        Mv = p.mv*obj.m_Pl;
        Mh = p.mh*obj.m_Pl;
        
        if p.n_recycle == 1
            phi0 = 0;
        else
            phi0_spacing = 8*Mh; % Separation between starting points on the same potential
            phi0 = phi0_spacing*(-floor(p.n_recycle/2):floor(p.n_recycle/2)-1+mod(p.n_recycle,2));
        end
        
        i_iter = 0;
        tic
        
        % Run n_iter iterations of the inflation simulation
        while true
            
            % Increment the counter at every loop
            if p.n_iter > 0, i_iter = i_iter + 1; end;
            if i_iter > abs(p.n_iter), break, end
            
            data_out = nan(1,17); % Initialize output array
            
            if (p.n_iter > 0 && mod(i_iter,1e4) == 0) || ...
               (p.n_iter < 0 && mod(i_iter,1e1) == 0)
                disp(num2str([i_iter toc]));
            end
            
            for goto = NaN % on break, goto record data
                
                record_flag = 0;
                
                %% Draw random potential and set initial conditions
                
                switch p.measure
                    
                    case 'A'
                        
                        % Generate dimensionless potential function f and
                        % compute potential values at starting points
                        if mod(i_iter,p.n_recycle) == 1
                            [ak,f{1:3}] = obj.gaussian_random_field_1D(p.kmax,p.gamma);
                            V0  = Mv^4*f{1}(phi0/Mh);
                            Vp0 = Mv^4*f{2}(phi0/Mh)/Mh;
                        end
                        
                        % Check if potential energy is already negative
                        ii = mod(i_iter-1,p.n_recycle)+1;
                        if V0(ii) < 0
                            data_out(2) = 1; break
                        else
                            data_out(2) = 0;
                        end
                        
                        % Look uphill for the (dim-less) potential peak
                        xpeak = find_phipeak(phi0(ii)/Mh,ak,1,0,1);
                        
                        % Start at maximum
                        % find_phipeak should ensure that the maximum is
                        % slightly toward the side we started on
                        if phi0(ii) < xpeak*Mh
                            phistart = xpeak*Mh - 1i; % To the left
                        else
                            phistart = xpeak*Mh + 1i; % To the right
                        end
                        Vstart  = [];
                        Vpstart = [];
                        
                    case 'B'
                        
                        % Generate dimensionless potential function f and
                        % compute potential values at starting points
                        if mod(i_iter,p.n_recycle) == 1
                            [ak,f{1:3}] = obj.gaussian_random_field_1D(p.kmax,p.gamma);
                            V0   = Mv^4*f{1}(phi0/Mh);
                            Vp0  = Mv^4*f{2}(phi0/Mh)/Mh;
                            Vpp0 = Mv^4*f{3}(phi0/Mh)/Mh^2;
                        end
                        
                        % Check if slow roll is valid at starting point
                        ii = mod(i_iter-1,p.n_recycle)+1;
                        if V0(ii) < 0
                            data_out(2) = 1;
                        elseif (Vp0(ii)/V0(ii)).^2/(16*pi/obj.m_Pl^2) > 1
                            data_out(2) = 2;
                        elseif abs(Vpp0(ii)/V0(ii))/(8*pi/obj.m_Pl^2) > 1
                            data_out(2) = 3;
                        else
                            data_out(2) = 0;
                        end
                        
                        % Abort if slow roll is invalid
                        if data_out(2) > 0, break, end
                        
                        phistart = phi0(ii);
                        Vstart   = V0(ii);
                        Vpstart  = Vp0(ii);
                        
                    case 'C'
                        
                        % Generate dimensionless potential function f and
                        % compute potential values at starting points
                        if mod(i_iter,p.n_recycle) == 1
                            [ak,f{1:3}] = obj.gaussian_random_field_1D(p.kmax,p.gamma);
                            V0   = Mv^4*f{1}(phi0/Mh);
                            Vp0  = Mv^4*f{2}(phi0/Mh)/Mh;
                            Vpp0 = Mv^4*f{3}(phi0/Mh)/Mh^2;
                        end
                        
                        % Check if slow roll is valid at starting point
                        ii = mod(i_iter-1,p.n_recycle)+1;
                        if V0(ii) < 0
                            data_out(2) = 1;
                        elseif (Vp0(ii)/V0(ii)).^2/(16*pi/obj.m_Pl^2) > 1
                            data_out(2) = 2;
                        elseif abs(Vpp0(ii)/V0(ii))/(8*pi/obj.m_Pl^2) > 1
                            data_out(2) = 3;
                        else
                            data_out(2) = 0;
                        end
                        
                        phistart = phi0(ii);
                        Vstart   = V0(ii);
                        Vpstart  = Vp0(ii);
                        
                end
                
                if obj.plotFlag
                    if isempty(obj.plotHandle)
                        obj.plotHandle = figure;
                    else
                        figure(obj.plotHandle)
                    end
                    x = linspace(-20,10,1001);
                    plot(x,f{1}(x)); hold on; plot(phi0/Mh,f{1}(phi0/Mh),'o'); hold off
                    pause(1)
                end
                
                %% Screen the potential
                
                rho_Lambda  = 0;
                
                % Find the first local minimum
                % Shift the potential up by rho_Lambda_thres so that
                % find_phistop doesn't give up if it goes slightly negative
                [xstop] = find_phistop(real(phistart)/Mh,ak,...
                    1, ... % Vscale
                    -p.rho_Lambda_thres^4/p.mv^4, ...
                    1, ... % phiscale
                    p.lambdascreenmode );
                
                near_minima = {};
                if p.rho_Lambda_thres > 0
                    
                    % Check if there is a minimum close to rho_Lambda_thres
                    
                    if isnan(xstop)
                        % Vacuum energy is too negative; abort
                        data_out(2) = 5;
                        break
                    end
                    
                    fstop = f{1}(xstop);
                    if abs(fstop) < p.rho_Lambda_thres^4/p.mv^4
                        % Vacuum energy is already close to a possible
                        % value for lambda in the starting basin
                        rho_Lambda = fstop*Mv^4;
                    else
                        
                        % Check if a neighboring vacuum is close to the
                        % observed value for lambda
                        
                        near_minima = {nan(1,p.n_tunnel_max),nan(1,p.n_tunnel_max)};
                        phipeak     = {nan(1,p.n_tunnel_max),nan(1,p.n_tunnel_max)};
                        phifv       = {xstop,xstop};
                        
                        found = false;
                        
                        % Randomize which direction is checked first
                        if randi(2) == 1
                            lr_order = [1 2];
                        else
                            lr_order = [2 1];
                        end
                        
                        rho_Lambda = zeros(1,2);
                        i_basin    = zeros(1,2);
                        
                        for lr = lr_order
                            for i_tunnel = 1:p.n_tunnel_max
                                [near_minima{lr}(i_tunnel),phipeak{lr}(i_tunnel)] = find_phinextmin(...
                                    phifv{lr}(i_tunnel),ak,1,0,1,1,2*(lr-1)-1,1);
                                if ...
                                        isnan(near_minima{lr}(i_tunnel)) || ...
                                        f{1}(near_minima{lr}(i_tunnel)) < 0 || ...
                                        diff(f{1}([near_minima{lr}(i_tunnel),phifv{lr}(i_tunnel)])) < 0
                                    break
                                elseif abs(f{1}(near_minima{lr}(i_tunnel))) < p.rho_Lambda_thres^4/p.mv^4
                                    rho_Lambda(lr) = f{1}(near_minima{lr}(i_tunnel))*Mv^4;
                                    i_basin(lr) = i_tunnel;
                                    found = true;    % Don't keep looking
                                    break
                                elseif i_tunnel < p.n_tunnel_max
                                    phifv{lr}(i_tunnel+1) = near_minima{lr}(i_tunnel);
                                end
                            end
                        end
                        
                        if found
                            if i_basin(1) == i_basin(2)
                                [~,ia] = sort(abs(rho_Lambda));
                            else
                                [~,ia] = sort(i_basin);
                            end
                            rho_Lambda = rho_Lambda(ia(1));
                        else
                            break % Skip to record data
                        end
                        
                    end
                    
                    if obj.plotFlag
                        if isempty(obj.plotHandle)
                            obj.plotHandle = figure;
                        else
                            figure(obj.plotHandle)
                        end
                        x = linspace(-8,8,1001)+real(phistart)/Mh;
                        plot(x,f{1}(x)); hold on; 
                        plot(real(phistart)/Mh,f{1}(real(phistart)/Mh),'o','LineWidth',1); 
                        plot(xstop,f{1}(xstop),'x','LineWidth',1); 
                        if ~isempty(near_minima) && ~all(isnan([near_minima{:}]))
                            plot([phipeak{:}],f{1}([phipeak{:}]),'^');
                            plot([near_minima{:}],f{1}([near_minima{:}]),'v','LineWidth',1);
                            pause(5)
                        end
                        hold off
                    end
                    
                elseif p.lambdascreenmode
                    
                    if isnan(xstop)
                        % Vacuum energy is negative; abort
                        data_out(2) = 5;
                        break
                    end
                    
                end
                
                %% Simulate inflation with tunneling events
                
                record_flag = 1;
                
                Vstop       = [];
                valid_basin = [];
                
                phi = nan(1,5);
                data_out(4) = false;
                
                goto_break = false;
                for i_tunnel = 0:p.n_tunnel_max
                    
                    %% Check for false-vacuum tunneling
                    
                    if i_tunnel > 0
                        
                        if isempty(Vstop), Vstop = V(phi(end)); end
                        
                        % Look for an instanton tunneling solution
                        [phitunnel,log_tunnel_rate,flag_hawking_moss,flag_fv_eternal] = ...
                            obj.check_false_vacuum_eternal(...
                            ak,f,phi,V,Vp,Vpp,rho_Lambda,p.n_tunnel_max+1-i_tunnel,Vstop,p.flag_screen_basin);
                        
                        if isnan(phitunnel) % No tunneling
                            if data_out(2) == 4 % No end of inflation
                                % Can't produce observables
                                goto_break = true;
                            end
                            break
                        end
                        
                        record_flag = 2;
                        
                        data_out(4) = data_out(4) || flag_fv_eternal;
                        % Record smallest tunneling rate
                        data_out(5) = min(data_out(5),log_tunnel_rate);
                        
                        if flag_hawking_moss
                            if phitunnel < phi(end)
                                % Start just to the left (-) of the maximum
                                phistart = phitunnel - 1i;
                            else
                                % Start just to the right (+) of the maximum
                                phistart = phitunnel + 1i;
                            end
                        else
                            phistart = phitunnel;
                        end
                        
                    end
                    
                    %% Simulate slowroll
                    % If tunneling occurs, simulate slow roll from the new
                    % starting point on other side of potential barrier
                    
                    it = 1 + i_tunnel;
                    [phi(it,:),status_sr(it,1),mv_sr(it,1),Ntotal(it,1),V,Vp,Vpp] = ...
                        obj.simulate_slowroll(f,phistart,rho_Lambda,Vstart,Vpstart);
                    
                    % Set phi(5) = phistop
                    if i_tunnel == 0
                        phi(it,5) = xstop*Mh;
                    else
                        [new_xstop] = find_phistop(real(phistart)/Mh,ak,...
                            1, ... % Vscale
                            rho_Lambda/Mv^4, ... % rho_Lambda
                            1, ... % phiscale
                            0 );   % lambdascreenmode
                        phi(it,5) = new_xstop*Mh;
                    end
                    
                    if V(phi(it,5)) <= 0
                        valid_basin(it) = false; %#ok<AGROW>
                        break
                    end
                    
                    if ~isnan(phi(it,4)) && ~(Ntotal(it) < p.Nafter)
                        valid_basin(1:it) = true;
                    else
                        valid_basin(it) = false; %#ok<AGROW>
                    end
                    
                end
                
                % Tunneling index of last valid basin
                i_last = find(valid_basin,1,'last');
                if isempty(i_last), i_last = length(mv_sr); end
                
                data_out(1) = mv_sr(i_last);
                data_out(2) = status_sr(i_last);
                data_out(3) = Ntotal(i_last);
                
                if goto_break
                    break % No end of inflation
                elseif ~(data_out(3) >= p.Nafter)
                    break % Not enough e-foldings
                end
                
                if obj.plotFlag
                    if isempty(obj.plotHandle)
                        obj.plotHandle = figure;
                    else
                        figure(obj.plotHandle)
                    end
                    x = linspace(-4*Mh,4*Mh,1001)+phi(i_tunnel+1,2);
                    plot(x,V(x)); hold on;
                    plot(phi(i_tunnel+1,1),V(phi(i_tunnel+1,1)),'^','LineWidth',1);
                    plot(phi(i_tunnel+1,2),V(phi(i_tunnel+1,2)),'o','LineWidth',1);
                    plot(phi(i_tunnel+1,3),V(phi(i_tunnel+1,3)),'p','LineWidth',1);
                    plot(phi(i_tunnel+1,4),V(phi(i_tunnel+1,4)),'*','LineWidth',1);
                    plot(phi(i_tunnel+1,5),V(phi(i_tunnel+1,5)),'x','LineWidth',1);
                    hold off
                    pause(3)
                end
                
                %% Check for stochastic and topological eternal inflation
                
                record_flag = 3;
                
                % Increment the counter only for successful inflation
                if p.n_iter < 0, i_iter = i_iter + 1; end
                
                for i_tunnel = 0:size(phi,1)-1
                    
                    it = 1+i_tunnel;
                    
                    % Set phi(1) = phipeak
                    if isreal(phistart)
                        phi(it,1) = find_phipeak(phi(it,2)/Mh,ak,1,rho_Lambda/Mv^4,1)*Mh;
                    else
                        phi(it,1) = phi(it,2);
                    end
                    
                    [numStochEpochs,NSinceStoch] = obj.check_stochastic_eternal(...
                        V,Vp,Vpp,phi(it,:));
                    flag_topological_eternal = obj.check_topological_eternal(...
                        V,Vp,Vpp,phistart,phi(it,1));
                    
                    data_out(6) = max(0,data_out(6)) + max(0,numStochEpochs);
                    data_out(7) = NSinceStoch; % Only keep the value from the last basin
                    data_out(8) = max(0,data_out(8)) + max(0,flag_topological_eternal);
                    
                end
                
                %% Compute observables
                
                % Construct V'''(phi)
                f{4} = function_handle.empty;
                [~,f{1:4}] = obj.gaussian_random_field_1D(p.kmax,p.gamma,ak);
                Vppp = build_potential(f,p.mv*obj.m_Pl,p.mh*obj.m_Pl,3);
                
                Nbefore = Ntotal - p.Nafter; % e-folds before horizon crossing
                observables = obj.compute_observables(V,Vp,Vpp,Vppp,phi(end,:),Nbefore,obj.m_Pl);
                
                data_out(9)  = observables.Q(1);
                data_out(10) = observables.r(1);
                data_out(11) = observables.n_s(1);
                data_out(12) = observables.alpha(1);
                data_out(13) = observables.n_t(1);
                data_out(14) = observables.dlgrho(1);
                data_out(15) = observables.lgOk(1);
                data_out(16) = observables.rho_Lambda(1)/obj.m_Pl^4; % Save dimensionless value
                
            end % for goto
            
            %% Record results
            
            % Record Flags
            % 1     Potential is valid; simulated slow roll
            % 2     Checked for false vacuum tunneling
            % 3     Found a valid basin with sufficient N_e
            
            % Idx   Output          Type        Record Flag
            % 1     mv              (float,4)   1,2,3
            % 2  	status          (int)       1,2,3
            % 3     Ntotal          (float,2)   1,2,3
            % 4  	flag_fv_eternal (bool)      2,3
            % 5  	log_tunnel_rate (float,4)   2,3
            % 6  	numStochEpochs  (int)       3
            % 7  	NSinceStoch     (float,2)   3
            % 8  	numTopolEpochs  (int)       3
            % 9  	Q               (float,4)   3
            % 10 	r               (float,4)   3
            % 11 	n_s             (float,4)   3
            % 12 	alpha           (float,4)   3
            % 13 	n_t             (float,4)   3
            % 14 	dlgrho          (float,4)   3
            % 15 	lgOk            (float,4)   3
            % 16 	rho_Lambda      (float,4)   3
            
            % Open output file for recording results
            % Record various amounts of data based on how far code ran
            
            switch record_flag
                case 1
                    fid = fopen(p.outfile,'at');
                    fprintf(fid,'%d,%.4G,%d,%.2G\r\n',...
                        [1 data_out(1:3)]);
                    fclose(fid);
                case 2
                    fid = fopen(p.outfile,'at');
                    fprintf(fid,'%d,%.4G,%d,%.2G,%d,%.4G\r\n',...
                        [2 data_out(1:5)]);
                    fclose(fid);
                case 3
                    fid = fopen(p.outfile,'at');
                    fprintf(fid,'%d,%.4G,%d,%.2G,%d,%.4G,%d,%.2G,%d,%.4G,%.4G,%.4G,%.4G,%.4G,%.4G,%.4G\r\n',...
                        [3 data_out(1:16)]);
                    fclose(fid);
            end
            
        end % for i_iter
        
    end
    
end

%% Slowroll

methods (Access = protected)
    
    function [phi,status,mv,Ntotal,V,Vp,Vpp] = simulate_slowroll(obj,f,phistart,rho_Lambda,Vstart,Vpstart)
        % Simulate slow roll inflation and collect data
        %
        % Inputs
        %   f   A cell array of function handles
        %           {V(phi), V'(phi), V''(phi), V'''(phi)}
        %
        % Outputs
        %   datastruct  A structure containing observables and diagnostics
        
        if nargin < 4, rho_Lambda = 0;  end
        if nargin < 5 || isempty(Vstart),  Vstart = [];     end
        if nargin < 6 || isempty(Vpstart), Vpstart = [];    end
        
        % If starting at a peak of the potential, 
        % the direction in which the inflaton will evolve
        fall_direction = sign(imag(phistart));
        phistart = real(phistart); % Initial value of the inflaton
        
        mv = obj.parameters.mv;
        mh = obj.parameters.mh;
        
        Ntotal = nan;
        
        phi = nan(1,5); % [phipeak, phistart, phiexit, phiend, phistop]
        phi(2) = phistart;
        
        kappa = 8*pi/obj.m_Pl^2;
        
        V0 = build_potential(f,mv*obj.m_Pl,mh*obj.m_Pl);
        V = @(x) V0(x) - rho_Lambda;
        Vp = []; Vpp = [];
        
        if nargin < 5, Vstart = V(phistart); end
        if nargin < 6
            Vp = build_potential(f,mv*obj.m_Pl,mh*obj.m_Pl,1);
            Vpstart = Vp(phistart);
        end
        
        %% Simulate slow roll
        
        for shift_potential = 0:uint8(obj.parameters.fixQ)
            
            [Vp,Vpp] = build_potential(f,mv*obj.m_Pl,mh*obj.m_Pl,[1 2]);
            
            %% Identify phiend, phiexit, and phistop
            
            % Find the end of slow roll inflation
            [phiend,status] = obj.find_phiend(phistart,V,Vp,Vpp,...
                mh*obj.m_Pl,obj.m_Pl,Vstart,Vpstart);
            if status == 4 % Slow roll valid until minimum; no phiend
                phi(5) = phiend;
                return
            else
                phi(4) = phiend;
            end
            
            % Compute total number of slow roll e-folds
            dlna_dphi = @(phi) (-kappa*V(phi)./Vp(phi));
            if fall_direction ~= 0
                Ntotal = Inf; % Started at a maximum
            else
                points = linspace(phistart,phiend,max(10,abs(phistart-phiend)/mv/obj.m_Pl*10));
                Ntotal_trapz = trapz(points,dlna_dphi(points));
                if  Ntotal_trapz > 0.7*obj.parameters.Nafter && ...
                        Ntotal_trapz < 1.3*obj.parameters.Nafter
                    Ntotal = integral(dlna_dphi,phistart,phiend);
                else
                    Ntotal = Ntotal_trapz;
                end
            end
            
            % Find the value of the field 55 e-folds
            % before the end of slow roll inflation
            if Ntotal >= obj.parameters.Nafter
                
                % Integrate back from (phiend,Nafter) to (phiexit,0)
                dphi_dlna = @(N,phi) (-Vp(phi)./V(phi)/kappa).';
                [~,phiexit] = ode23(dphi_dlna,[obj.parameters.Nafter 1 0],phiend);
                phi(3) = phiexit(3);
                
            end
            
        end
        
    end
    
end

%% Checks for eternal inflation

methods (Access = protected)
    
    function [phitunnel,log_tunnel_rate,flag_hawking_moss,flag_fv_eternal] = check_false_vacuum_eternal(obj,ak,f,...
            phi,V,Vp,Vpp,rho_Lambda,n_tunnel_remaining,Vstop,flag_screen_basin)
        % Handle false vacuum tunneling.
        
        if nargin < 7 || isempty(rho_Lambda), rho_Lambda = 0; end
        if nargin < 8, n_tunnel_remaining = 1; end
        if nargin < 9, Vstop = V(phi(5)); end
        if nargin < 10, flag_screen_basin = true; end
        
        flag_rescale_potential = false;
        
        mv = obj.parameters.mv;
        mh = obj.parameters.mh;
        
        kappa = 8*pi/obj.m_Pl^2;
        hbar  = 1;
        
        phiscale = mh*obj.m_Pl;
        
        flag_fv_eternal = true;
        log_stable_rate_cutoff = -1e2;
        log_tunnel_rate = nan;
        
        %% Locate nearest minima
        
        phistop = phi(5); % Location of false vacuum
        
        near_minima = {NaN,NaN};
        
        % Find minima around phistop
        Vpstop = Vp(phistop);
        
        phifv = {phistop,phistop};
        phipeak = {nan,nan};
        Vfv   = {Vstop,Vstop};
        Vpfv  = {Vpstop,Vpstop};
        
        %% Check for tunneling to nearby minima
        
        new_phistart = nan(1,2);
        flag_hawking_moss = false(1,2);
        B = inf(1,2);
        
        for lr = 1:2 % Left and right neighbor basins
            
            if true %flag_screen_basin
                % Look ahead up to next n_tunnel_remaining adjacent basins. If
                % we can't get enough e-foldings on the other side the
                % potential barrier, don't bother computing the tunneling rate
                for i_tunnel = 1:n_tunnel_remaining
                    
                    lr_continue = true;
                    
                    [near_minima{lr}(i_tunnel),phipeak{lr}(i_tunnel)] = find_phinextmin(...
                        phifv{lr}(i_tunnel),ak,mv^4*obj.m_Pl^4,rho_Lambda,phiscale,1,2*(lr-1)-1,1);
                    
                    if isnan(near_minima{lr}(i_tunnel)) || ...
                            V(near_minima{lr}(i_tunnel)) > Vfv{lr}(i_tunnel)
                        break
                    end
                    
                    if ~flag_screen_basin
                        lr_continue = false;
                        break;
                    end
                    
                    %% Find barrier edge location
                    % Use this starting point to compute the maximum amount of
                    % inflation that could occur in the next potential basin
                    
                    phi_tol = abs(phifv{lr}(i_tunnel)-near_minima{lr}(i_tunnel))*1e-10;
                    
                    phimin = phifv{lr}(i_tunnel);
                    phimax = near_minima{lr}(i_tunnel);
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
                            fun_vals = Vfv{lr}(i_tunnel)-V(phi_range);
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
                    
                    Vbar = V(phibar);
                    if (Vp(phibar)/Vbar).^2/(2*kappa) > 1
                        status_bar = 1;
                    elseif abs(Vpp(phibar)/Vbar)/(kappa) > 1
                        status_bar = 2;
                    else
                        status_bar = 0;
                    end
                    
                    %% Find the start of inflation in the neighboring basin
                    
                    if status_bar == 0 % Inflation at barrier edge
                        phistart = phibar;
                    else
                        % Look for inflation further down the slope
                        [phistart,status_start] = obj.find_phistart_downhill(...
                            phibar,V,Vp,Vpp,phiscale,obj.m_Pl);
                    end
                    
                    %% Check if there is enough inflation
                    
                    if status_bar == 0 || status_start == 0 % Found inflation
                        
                        % Find the end of slow roll inflation in neighbor basin
                        [phiend,status_end] = obj.find_phiend(phistart,V,Vp,Vpp,...
                            phiscale,obj.m_Pl,[],[],false,(i_tunnel < n_tunnel_remaining));
                        
                        Ntotal = nan;
                        
                        if status_end == 4
                            % End of inflation does not occur before reaching the
                            % next local minimum - cannot produce an observable
                            % universe in that basin of the potential
                            
                            if i_tunnel == n_tunnel_remaining
                                break
                            end
                            
                        else
                            
                            % Compute total number of slow roll e-folds
                            dlna_dphi = @(phi) (kappa*V(phi)./Vp(phi));
                            points = linspace(phistart,phiend,max(10,abs(phistart-phiend)/mh/obj.m_Pl*10));
                            Ntotal_trapz = trapz(points,dlna_dphi(points));
                            if  Ntotal_trapz > 0.7*obj.parameters.Nafter
                                Ntotal = integral(dlna_dphi,phistart,phiend);
                            else
                                Ntotal = Ntotal_trapz;
                            end
                            
                        end
                        
                    else
                        
                        status_end = -1;
                        Ntotal = 0;
                        
                    end
                    
                    if status_end ~= 4
                        % Reheating occurs
                        % A candidate observable universe
                        
                        if Ntotal >= obj.parameters.Nafter
                            
                            % Tunneling basin is viable
                            lr_continue = false;
                            break
                            
                        else
                            
                            %% Pre-compute Hawking-Moss tunneling rate
                            %  Check if a Hawking-Moss instanton transition can
                            %  give us inflation near the local maximum
                            
                            w_top = abs(kappa/3*V(phipeak{lr}(i_tunnel)))^0.5;
                            R = pi/w_top/2; Y = [phipeak{lr}(i_tunnel),0,1/w_top,0,-w_top];
                            B_HM = FalseVacuumInstanton.find_tunneling_suppression_static(...
                                3,kappa,V,phifv{lr}(i_tunnel),R,Y);
                            
                            if -B_HM >= log_stable_rate_cutoff
                                % Hawking-Moss instanton is viable
                                lr_continue = false;
                                break
                            elseif i_tunnel == n_tunnel_remaining
                                % No more chances to tunnel
                                break
                            end
                            
                        end
                        
                    end
                    
                    if i_tunnel < n_tunnel_remaining
                        phifv{lr}(i_tunnel+1) = near_minima{lr}(i_tunnel);
                        Vfv{lr}(i_tunnel+1)   = V(phifv{lr}(i_tunnel+1));
                        Vpfv{lr}(i_tunnel+1)  = Vp(phifv{lr}(i_tunnel+1));
                    end
                    
                    lr_continue = false;
                    
                end
                
                if lr_continue, continue, end
                
                disp('Computing instanton...');
                
            end
            
            %% Do full instanton calculation
            
            xtol         = 1e-4;
            phitol       = 1e-4;
            thinCutoff   = 1e-2;
            B_cutoff     = 1e4;
            try
            if flag_rescale_potential
                
                % Define potential with rescaled mv = 1
                [V_rescale,Vp_rescale,Vpp_rescale] = build_potential(f,1,mh*obj.m_Pl);
                f_Lambda = rho_Lambda/mv^4/obj.m_Pl^4;
                V_rescale = @(x) V_rescale(x) - f_Lambda;
                
                % Initialize instanton solver
                fvi = FalseVacuumInstanton(...
                    'V',            V_rescale,...
                    'dV',           Vp_rescale,...
                    'd2V',          Vpp_rescale,...
                    'm_Pl',         obj.m_Pl,...
                    'phi_metaMin',  phistop,...
                    'phi_absMin',   near_minima{lr}(1) );
                
            else
                
                % Initialize instanton solver
                fvi = FalseVacuumInstanton(...
                    'V',            V,...
                    'dV',           Vp,...
                    'd2V',          Vpp,...
                    'm_Pl',         obj.m_Pl,...
                    'phi_metaMin',  phistop,...
                    'phi_absMin',   near_minima{lr}(1) );
                
            end
            catch me
                continue
            end

            try % Solve for instanton profile
                [R,Y,~] = fvi.find_profile([],xtol,phitol,thinCutoff);
            catch me
                switch me.identifier
                    case 'FalseVacuumInstanton:StableFalseVacuum'
                        continue % No tunneling
                    case 'FalseVacuumInstanton:IntegralDiverged'
                        continue % Integration failed; assume no tunneling
                    otherwise
                        continue; %rethrow(me);
                end
            end
            
            if isscalar(R), flag_hawking_moss(lr) = true; end
            
            if isnan(Y(1,1))
                continue
            else
                new_phistart(lr) = Y(1,1); % Field value at center of bubble
            end
            
            % Get tunneling suppression rate B = -log(\lambda)
            % using appropriately scaled radial coordinate
            % B = S_bubble - S_background
            if flag_rescale_potential
                B(lr) = fvi.find_tunneling_suppression(R,Y)/mv^4/obj.m_Pl^4;
            else
                B(lr) = fvi.find_tunneling_suppression(R,Y);
            end
            
        end
        
        %% Choose bubble with the larger tunneling rate
        
        if all(isnan(new_phistart))
            phitunnel = NaN;
            return
        end
        
        % B has units of kg*m
        [log_tunnel_rate,imax] = max(-B/hbar);
        
        %% Find the new value of phi after tunneling
        
        if log_tunnel_rate > log_stable_rate_cutoff
            % Tunneling is fast enough
            phitunnel = new_phistart(imax);
            flag_hawking_moss = flag_hawking_moss(imax);
        else
            % Tunneling is too slow
            phitunnel = NaN;
            flag_hawking_moss = false;
        end
        
        %% Determine if inflation is eternal
        
%         H4 = (kappa/3*V(phistop))^2;
        if log_tunnel_rate < log(9/4/pi);
            % Tunneling is slow enough for eternal inflation
            flag_fv_eternal = true;
        else
            % Tunneling is too fast for eternal inflation
            flag_fv_eternal = false;
        end
        
    end
    
    function [numStochEpochs,NSinceStoch] = check_stochastic_eternal(obj,V,Vp,~,phi)
        
        p = obj.parameters;
        
        numStochEpochs = nan;
        NSinceStoch = nan;
        
        phipeak  = phi(1);
        phistart = phi(2);
        phiexit  = phi(3);
        phiend   = phi(4);
        
        kappa = 8*pi/obj.m_Pl^2;
        hbar = 1;
        
        phiscale = p.mh*obj.m_Pl;
        
        %% find phistoch
        
        Vstart  = V(phistart);
        Vpstart = Vp(phistart);
        
        sgn    = sign(Vp(phistart));
        phimin = phiscale^2 * abs(Vpstart./Vstart);
        dphi   = -0.001*sgn*max(0.01*phiscale,min(phimin,phiscale));
        
        % SEIC = Stochastic Eternal Inflation Criterion
        
        phibreak = nan; % Field values where SEIC changes sign
        off2on = nan; % True if SEIC goes from (+) to (-) moving downhill
        phibreak_last = nan;
        
        % Positive if stochastic eternal
        sgn_sei = sign(hbar^0.5*(kappa*Vstart)^1.5 - 2*pi*sqrt(3)*(0.607)*abs(Vpstart));
        ind_sei = 1;
        
        step = 10^(1/16);
        ind = 1; ii = 1;
        batch = 40;
        phi = phistart*ones(1,batch);
        while sgn*(phi(ii)-phiend) > 0
            if mod(ind,batch) == 1
                phi_last = phi(end);
                phi = phi_last + cumsum(dphi*step.^(ind-1:ind+batch-2));
                Vend   = V(phi);
                Vpend  = Vp(phi);
            end
            ii = mod(ind-1,batch)+1;
            if sgn_sei*(hbar^0.5*(kappa*Vend(ii))^1.5 - 2*pi*sqrt(3)*(0.607)*abs(Vpend(ii))) < 0
                if ii > 1
                    phibreak_last(ind_sei) = phi(ii-1);
                else
                    phibreak_last(ind_sei) = phi_last;
                end
                phibreak(ind_sei) = phi(ii);
                if sgn_sei == -1
                    off2on(ind_sei) = 1;
                else
                    off2on(ind_sei) = 0;
                end
                sgn_sei = -sgn_sei;     % Look for next change in status
                phi = phi(ii); ii = 1;  % Reset the rate of advance
                ind_sei = ind_sei + 1;
                ind = 0;
            end
            ind = ind + 1;
        end
        
        % Reached phiend without satisfying eternal inflation criterion
        if isnan(phibreak(1)), return, end
        
        % Find field values at breakdown precisely using binary search
        fun = @(x) (hbar^0.5*kappa*V(x)).^1.5 - 2*pi*sqrt(3)*(0.607)*abs(Vp(x));
        for i_break = 1:length(phibreak)
            
            % Initial bounds
            phimax = phibreak(i_break);
            phimin = phibreak_last(i_break);
            phisep = (phimax-phimin);
            
            phibreak(i_break) = 0.5*(phimin+phimax);
            
            nbits = 5;
            ind = 1;
            while abs(phisep)*2^(1-ind) > abs(dphi)
                
                % Precompute a binary tree of function values
                if mod(ind,nbits-1) == 1
                    phisep = (phimax-phimin);
                    phi = (phimin + phisep*(2^-nbits)):(phisep*(2^-nbits)):(phimax - phisep*(2^-nbits));
                    fun_vals = fun(phi);
                    ind = 1;
                    ii = 2^(nbits-1);
                end
                
                % Update bounds
                if fun_vals(ii) > 0
                    phimax = phi(ii);
                    ii = ii - 2^(nbits-1-ind);
                else
                    phimin = phi(ii);
                    ii = ii + 2^(nbits-1-ind);
                end
                
                % Update value
                phibreak(i_break) = 0.5*(phimin+phimax);
                ind = ind + 1;
                
            end
        
        end
        
        if off2on(1) == 0 % SEI is satisfied at phistart
            % Require that SEI ends when inflation ends, if not sooner
            phibreak = [phistart phibreak];
            off2on   = [1        off2on];
        end
        
        if off2on(end) == 1 % SEI is satisfied at phiend
            % Require that SEI ends when inflation ends, if not sooner
            phibreak(end+1) = phiend;
            off2on(end+1)   = 0;
        end
        
        dlna_dphi = @(phi) -V(phi)./Vp(phi); % Integrate this to find N_e
        
        % Check that at least one e-fold elapses with SEI valid
        for i = find(off2on,1)
            N = integral(@(x) dlna_dphi(x),phibreak(i),phibreak(i+1));
            if abs(N) < 1
                % SEI epoch is too short; don't count it
                phibreak(i:i+1) = [];
                off2on(i:i+1) = [];
            end
        end
        
        if isempty(phibreak), return, end
        
        numStochEpochs = nnz(~off2on);
        
        % Compute # of e-folds after past SEI breakdown and before phiexit
        if isnan(phiexit)
            return % not applicable
        elseif (hbar^0.5*kappa*V(phiexit)).^1.5 - 2*pi*sqrt(3)*(0.607)*abs(Vp(phiexit)) > 0
            NSinceStoch = 0; % Eternal at phiexit
        elseif any((phibreak(~off2on)-phiexit)*Vp(phiexit) > 0)
            % Closest breakdown point of SEI to phiexit
            phibreak_last = phibreak((phibreak-phiexit)*Vp(phiexit) > 0 & ~off2on);
            [~,imin] = min(abs(phibreak_last-phiexit));
            % Integrate e-foldings
            NSinceStoch = integral(@(x) dlna_dphi(x),phibreak_last(imin),phiexit);
        end
        
        return
        
        %% Stochastic Eternal inflation
        
% %         % Here \kappa = 1
% %         stochasticEIC = @(x) (kappa*V(x)).^(3/2) > 2*pi*sqrt(3)*(0.607)*abs(Vp(x)); % stochasticEIC > 0 -> Eternal
% %         
% %         % Find the value of the field 55 e-folds
% %         % before the end of slow roll inflation
% %         
% %         lna_steps = 0:floor(Ntotal);
% %         if length(lna_steps) == 1, return, end
% %         [~,phi_N] = ode45(dphi_dlna,lna_steps,phistart);                 % Phi values separated by 1 e-fold of inflation
% %         phi_N = sort([phi_N; phiend; linspace(phistart,phimax,20).']).'; % Add search points between phistart and phimax
% %         if phimax > phistart, phi_N = fliplr(phi_N); end                 % Order [phimax; ... phistart; ... phiend]
% %         phi_N(isnan(phi_N)) = [];
% %         
% %         % Find where stochastic eternal inflation starts/ends
% %         % Should always be eternal at phimax
% %         intv_ind = find(diff(stochasticEIC(phi_N)));
% %         off2on_sei   = feval(@(x) x(intv_ind) > 0, diff(stochasticEIC(phi_N)));
% %         phibreak_sei = arrayfun(@(i) fzero(stochasticEIC,phi_N(i:i+1)),intv_ind);
% %         
% %         if isempty(phibreak_sei)
% %             phibreak_sei = phimax;
% %             off2on_sei = 0;
% %         end
        
% %         % Second derivative check for stochastic eternal inflation near the
% %         % maximum. Assume potential is locally mirror symmetric about maximum,
% %         % or we condition on the field falling only toward phistart.
% %         phidev = phimax-phibreak_sei(1); Hmax = sqrt(V(phimax)/3);
% %         if abs(erf(abs(phidev/(Hmax/2/pi))/sqrt(2))) > exp(-3)
% %             numStochEpochs = uint8(1 + nnz(off2on_sei));
% %         else
% %             stochasticEIC = @(x) stochasticEIC(x) & ...
% %                 ~(min(phimax,phibreak_sei(1)) <= x && x <= max(phimax,phibreak_sei(1)));
% %             phibreak_sei(1) = []; off2on_sei(1) = [];
% %             numStochEpochs = uint8(nnz(off2on_sei));
% %         end
        
% %         % Compute # of e-folds after past SEI breakdown and before phiexit
% %         if stochasticEIC(phiexit)
% %             NSinceStoch = 0; % Eternal at phiexit
% %         elseif any((phibreak_sei(~off2on_sei)-phiexit)*Vp(phiexit) > 0)
% %             phib = phibreak_sei((phibreak_sei-phiexit)*Vp(phiexit) > 0 & ~off2on_sei);
% %             [~,imin] = min(abs(phib-phiexit));
% %             NSinceStoch = integral(@(x) dlna_dphi(x).*slowroll(x),...
% %                 phib(imin),phiexit);
% %         else
% %             NSinceStoch = NaN; % Never eternal before phiexit
% %         end
        
        % Compute # of e-folds eternal after phistart and before phiend
%         if strcmpi(p.measure,'a') && stochasticEIC(phistart)
%             NStochastic = Inf;
%         else
%             NStochastic = nan;
% %             datastruct.NStochastic = integral(@(x) ...
% %                 dlna_dphi(x).*stochasticEIC(x).*slowroll(x),phistart,phiend);
%         end
        
    end
    
    function [flag_topological_eternal] = check_topological_eternal(obj,V,Vp,Vpp,phistart,phipeak)
        
        p = obj.parameters;
        
        flag_topological_eternal = nan;
        
        %% Topological Eternal Inflation
        
        kappa = 8*pi/obj.m_Pl^2;
        hbar  = 1;
        
        % Determine whether quantum fluctuations could result in at least
        % one Hubble volume descending toward a different minimum of V(phi)
        if isreal(phistart) % Not Measure A
            H = sqrt(kappa/3*V(phistart));
            dphi = hbar^0.5*H/2/pi; % Amplitude of quantum fluctuation: P_\phi(k) = (dphi)^2 = \hbar * (H/2/\pi)^2
            Dphi = -Vp(phistart)/(3*H^2); % Classical field excursion
            if 1/2*erfc(abs(phistart+Dphi-phipeak)/(sqrt(2)*dphi)) < exp(-3)
                % Not close enough to maximum
                return
            end
        end
        
        % Assuming a domain wall forms around maximum, does it persist?
        
        % Find value of phi at the domain wall boundary 
        % where inflation breaks down
        [phiedge,status] = obj.find_phiend(phipeak,V,Vp,Vpp,p.mh*obj.m_Pl,obj.m_Pl,[],[],false,false);
        if status == 1
            flag_topological_eternal = false;
            return
        elseif status == 4
            % Inflation doesn't end
            flag_topological_eternal = true;
            return
        end
        
        % Find value of phi that will descend to phiedge in time < t_H
        dphi_dlna = @(lna,phi) -Vp(phi)./V(phi)/kappa;
        [~,phiout] = ode23(dphi_dlna,[1,0.5,0],phiedge);
        phistar = phiout(3);
        
        if isnan(phistar), return, end
        
        % If the expansion of the inner portion of the domain wall
        % where (phipeak < phi < phistar) replaces loss of the outer
        % wall where (phistar < phi < phiedge), then inflation is eternal
        % (Problem is effectively 1D - the width of domain wall.)
        if abs(phistar-phipeak) > abs(phiedge-phipeak)*exp(-1)
            flag_topological_eternal = true;
        end
        
% %         phiedge_eps = fzero(@(phi) Vp(phi)./V(phi)/2/kappa - 1,phipeak);
% %         phiedge_eta = fzero(@(phi) abs(Vpp(phi)./V(phi)/kappa) - 1,phipeak);
% %         [~,icloser] = min(abs(phipeak-[phiedge_eps,phiedge_eta]));
% %         phiedge = feval(@(x) x(icloser),[phiedge_eps,phiedge_eta]);
        
    end
    
end

%% Efficient search functions

methods (Static)
    
    function [phiend,status] = find_phiend(phistart,V,Vp,Vpp,phiscale,m_Pl,Vstart,Vpstart,lambdascreenmode,precisephistop)
        % Find the value of phi at the end of slow roll inflation
        
        if nargin < 7 || isempty(Vstart), Vstart = V(phistart);     end
        if nargin < 8 || isempty(Vpstart), Vpstart = Vp(phistart);   end
        if nargin < 9, lambdascreenmode = false; end
        if nargin < 10, precisephistop = true; end
        
        kappa = 8*pi/m_Pl^2;
        
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
            elseif (Vpend(ii)/Vend(ii)).^2/(2*kappa) > 1
                status = 2;
            elseif abs(Vppend(ii)/Vend(ii))/kappa > 1
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
            case 1, fun = V;                                    % V = 0
            case 2, fun = @(x) (Vp(x)./V(x)).^2/(2*kappa) - 1;  % eps = 1
            case 3, fun = @(x) abs(Vpp(x)./V(x))/kappa - 1;     % eta = 1
            case 4, fun = Vp;                                   % Vp = 0
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
    
    function [phistart,status] = find_phistart_downhill(phiinit,V,Vp,Vpp,phiscale,m_Pl)
        % Find the value of phi at the end of slow roll inflation
        
        sgn    = sign(Vp(phiinit));
        phimin = phiscale^2 * abs(Vp(phiinit)./V(phiinit));
        dphi   = -0.001*sgn*max(0.01*phiscale,min(phimin,phiscale));
        
        kappa = 8*pi/m_Pl^2;
        
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
            if (Vpend(ii)/Vend(ii)).^2/(2*kappa) > 1
                status = 2;
            elseif abs(Vppend(ii)/Vend(ii))/(kappa) > 1
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
            case 2, fun = @(x) (Vp(x)./V(x)).^2/(2*kappa) - 1;  % eps = 1
            case 3, fun = @(x) abs(Vpp(x)./V(x))/kappa - 1;   % eta = 1
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
    
end

methods (Static)
    
    function [status] = compute_status(V,Vp,Vpp,phi,m_Pl)
        if V(phi) < 0
            status = 1;
        elseif (Vp(phi)/V(phi)).^2/(16*pi/m_Pl^2) > 1
            status = 2;
        elseif abs(Vpp(phi)/V(phi))/(8*pi/m_Pl^2) > 1
            status = 3;
        else
            status = 0;
        end
    end
    
    function [observables] = compute_observables(V,Vp,Vpp,Vppp,phi,Nbefore,m_Pl)
        % Inputs
        %   potential   Cell array of function handles
        %                   {V, V', V'', V'''}
        %   phi         Monotonic array of field values
        %                   [phistart, phiexit, phiend, phistop]
        %   Nbefore     Number of e-foldings between phistart and phiexit
        %
        % Outputs
        %   observables Structure of observables
        
        phiexit  = phi(3);
        V_exit = V(phiexit);
        
        kappa = 8*pi/m_Pl^2;
        
        % Dimensionless slow roll parameters at horizon exit scale
        eps = (Vp(phiexit)./V_exit).^2/(2*kappa);
        eta = Vpp(phiexit)./V_exit/kappa;
        xi2 = Vp(phiexit).*Vppp(phiexit)./V_exit^2/kappa^2;
        
        observables.Q           = sqrt(V_exit/(150*eps)*kappa*m_Pl^2)/pi; % dimensionless
        observables.r           = 16*eps;
        observables.n_s         = 1-6*eps+2*eta;
        observables.alpha       = 16*eps*eta - 24*eps^2 - 2*xi2;
        observables.n_t         = -2*eps;
        observables.dlgrho      = log10(V_exit/V(phi(4)));
        observables.lgOk        = log10(V(phi(2))/V_exit) - Nbefore*2/log(10);
        observables.rho_Lambda  = V(phi(5));
        
    end
    
end

% Deprecated
methods (Static)
    
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
    
    function [phipeak,Vppeak] = find_phipeak(phi0,V,Vp,phiscale,V0,Vp0)
        % Find the value of phi at the previous local maximum
        
        if nargin < 6, V0  = V(phi0);  end
        if nargin < 7, Vp0 = Vp(phi0); end
        
        sgn    = sign(Vp0);
        phimin = phiscale^2 * abs(Vp0./V0);
        dphi   = -0.001*sgn*max(0.01*phiscale,min(phimin,phiscale));
        
        % Take steps until passed local max 
        step = 1e2;
        ind = 1;
        batch = 30;
        phi = phi0*ones(1,batch);
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
        sgn_Vp  = sgn*Vp(phipeak);
        while abs(phimax-phimin) > abs(dphi) || sgn_Vp < 0
            if sgn_Vp > 0
                phimax = phipeak;
            else
                phimin = phipeak;
            end
            phipeak = 0.5*(phimin + phimax);
            Vppeak = Vp(phipeak);
            sgn_Vp = sgn*Vppeak;
        end
        
    end
    
    function [phinextmin,phipeak] = find_next_minimum(phistop,V,Vp,Vpp,phiscale,search_direction,Vstop,Vpstop,alsopeak)
        % Find the value of phi at the next local minimum
        
        if nargin < 6, search_direction = 1; end
        if nargin < 9, alsopeak = true; end
        
        if nargin < 7 || isempty(Vstop),  Vstop = V(phistop);   end
        if nargin < 8 || isempty(Vpstop), Vpstop = Vp(phistop); end
        
        V0 = Vstop;
        
        phimin = phiscale^2 * abs(Vpstop./Vstop);
        dphi   = 0.001*search_direction*max(0.01*phiscale,min(phimin,phiscale));
        
%         x = linspace(-5*phiscale,5*phiscale,1001);
        
        phipeak = nan;
        
        % Take steps until passed local max
        flag_crossed_maximum = false;
%         step = 3e3;
        step = 10^(1/16);
        ind = 1;
        batch = 25;
        phinextmin = phistop*ones(1,batch);
%         dphi*step.^(ind-1:ind+batch-2)
        while true
            if mod(ind,batch) == 1
                phinextmin_last = phinextmin(end);
                phinextmin = phinextmin_last + cumsum(dphi*step.^(ind-1:ind+batch-2));
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
                if alsopeak
                    phipeakmin = phinextmin(ii);
                    if ii > 1
                        phipeakmax = phinextmin(ii-1);
                    else
                        phipeakmax = phinextmin_last;
                    end
                end
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
        
        if alsopeak && V(phinextmin) < V0
            
            phimin = phipeakmin;
            phimax = phipeakmax;
            if phimin > phimax
                phitemp = phimin;
                phimin = phimax;
                phimax = phitemp;
            end
            phipeak = 0.5*(phimin + phimax);
            phisep = abs(phimax-phimin);
            
            nbits = 5;
            ind = 1;
            while abs(phisep)*2^(1-ind) > abs(dphi)
                if mod(ind,nbits-1) == 1
                    phisep = (phimax-phimin);
                    phi = (phimin + phisep*(2^-nbits)):(phisep*(2^-nbits)):(phimax - phisep*(2^-nbits));
                    fun_vals = -Vp(phi);
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
                phipeak = 0.5*(phimin+phimax);
                ind = ind + 1;
            end
            
        end
        
    end
    
    function [phistart,status] = find_phistart_uphill(phiinit,V,Vp,Vpp,phiscale,m_Pl)
        % Find the value of phi at the end of slow roll inflation
        
        sgn    = sign(Vp(phiinit));
        phimin = phiscale^2 * abs(Vp(phiinit)./V(phiinit));
        dphi   = 0.001*sgn*max(0.01*phiscale,min(phimin,phiscale));
        
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
            if (Vpend(ii)/Vend(ii)).^2/(16*pi/m_Pl^2) > 1
                status = 2;
            elseif abs(Vppend(ii)/Vend(ii))/(8*pi/m_Pl^2) > 1
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
    
end

%% Utilities

methods (Static)
    
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
    
    function [datastruct] = read_output_file(outfile,record_flags)
        
        fid = fopen(outfile,'r');
        
        %% Collect meta-data
        
        meta_line = fgets(fid); is = 1;
        
        [n_iter,~,~,is1] = sscanf(meta_line(is:end),'%E,',1); is = is+is1-1;
        [mv_0,~,~,is1]   = sscanf(meta_line(is:end),'%G,',1); is = is+is1-1;
        [mh,~,~,is1]     = sscanf(meta_line(is:end),'%G,',1); is = is+is1-1;
        [m_Pl,~,~,is1]    = sscanf(meta_line(is:end),'%G,',1); is = is+is1-1;
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
        
        %% Read output from simulations
        
        data = nan(round(abs(n_iter)/10),16);
        mbytes = 0;
        
        for i_ln = 1:abs(n_iter)
            
            % Query record flag for this line
            [record_flag] = fscanf(fid,'%d,',1);
            
            % Reached end of file?
            if isempty(record_flag), break, end
            
            bytes = ftell(fid);
            
            if bytes/1024000.0 > mbytes+1
                mbytes = floor(bytes/1024000.0);
                disp([num2str(mbytes) ' MB']);
            end
            
            % Only read data corresponding to queried record flags
            if nargin >= 2 && ~ismember(record_flag,record_flags)
                fgets(fid); continue % Skip this line
            end
            
            % Read a line of data
            switch record_flag
                case 1
                    data(i_ln,1:3)  = fscanf(fid,'%G,%d,%G\r\n',3);
                case 2
                    data(i_ln,1:5)  = fscanf(fid,'%G,%d,%G,%d,%G\r\n',5);
                case 3
                    data(i_ln,1:16) = fscanf(fid,'%G,%d,%G,%d,%G,%d,%G,%d,%G,%G,%G,%G,%G,%G,%G\r\n',16);
            end
            
        end
        
        data(isnan(data(:,1)),:) = [];
        
        datastruct = struct(...
            'n_iter',           n_iter,...
            'mv_0',             mv_0,...
            'mh',               mh,...
            'm_Pl',             m_Pl,...
            'kmax',             kmax,...
            'gamma',            gamma,...
            'measure',          measure,...
            'n_tunnel_max',     n_tunnel_max,...
            'lambdascreen',     logical(lambdascreen),...
            'rho_Lambda_thres', rho_Lambda_thres,...
            'fixQ',             logical(fixQ),...
            'Nafter',           Nafter,...
            'seed',             seed,...
            'n_recycle',        n_recycle,...
            'data',             data);
        
        fclose(fid);
        
    end
    
    function plot_histograms(data,savename)
        
        if nargin < 2, savename = ''; end
        save_fig = [savename '_%s'];
        
        observables.numStochEpochs = data(:,6);
        observables.NSinceStoch = data(:,7);
        
        observables.Q       = data(:,9);
        observables.n_s     = data(:,11);
        observables.alpha   = data(:,12);
        observables.n_t     = data(:,13);
        observables.lgOk    = data(:,15);
        
%         save(sprintf(save_fig,'observables'),'observables');
        
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
%         close(gcf)
        
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
%         close(gcf)
        
        % alpha - Running of the spectral index
        fieldname = 'alpha';
        figure, hist([observables.(fieldname)],-0.51:0.02:0.51,hist_color)
        title('Running: \alpha')
%         set(gca,'XLim',[-0.5,0.5]);
        h = findobj(gca,'Type','patch');
        h.FaceColor = hist_color; boldify
        if ~isempty(savename)
            hgexport(gcf,'-clipboard',style,'applystyle',true); drawnow
            savefig(sprintf(save_fig,fieldname));
            set(gcf,'PaperPositionMode','auto');
            print(sprintf(save_fig,fieldname),'-dpng','-r0');
        end
%         close(gcf)
        
        % n_t - Tensor spectral index
        fieldname = 'n_t';
        figure, hist([observables.(fieldname)],-0.051:0.001:0.01,hist_color)
        title('Tensor Spectral Index: n_t')
%         set(gca,'XLim',[-0.05,0]);
        h = findobj(gca,'Type','patch');
        h.FaceColor = hist_color; boldify
        if ~isempty(savename)
            hgexport(gcf,'-clipboard',style,'applystyle',true); drawnow
            savefig(sprintf(save_fig,fieldname));
            set(gcf,'PaperPositionMode','auto');
            print(sprintf(save_fig,fieldname),'-dpng','-r0');
        end
%         close(gcf)
        
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
%         close(gcf)
        
        % Number of e-foldings since eternal inflation
        fieldname = 'NSinceStoch';
        val = log10([observables.(fieldname)]);
        figure, hist(val(isfinite(val)),12,hist_color)
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
%         close(gcf)
        
        % Number of eternal inflation epochs
        fieldname = 'numStochEpochs';
        val = [observables.(fieldname)];
        figure, hist(val(val>0),10,hist_color)
        title('Number of eternal inflation epochs')
        h = findobj(gca,'Type','patch');
        h.FaceColor = hist_color; boldify
        set(gca,'XTick',[1 2 3])
        if ~isempty(savename)
            hgexport(gcf,'-clipboard',style,'applystyle',true); drawnow
            savefig(sprintf(save_fig,fieldname));
            set(gcf,'PaperPositionMode','auto');
            print(sprintf(save_fig,fieldname),'-dpng','-r0');
        end
%         close(gcf)
        
    end
    
end

%% Constructor

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
                            isfinite(val.n_iter) && val.n_iter ~= 0 && ...
                            mod(val.n_iter,1) == 0 && isreal(val.n_iter)
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.n_iter must be a real, finite, positive integer.');
                    end
                    
                case 'mv'
                    
                    if isnumeric(val.mv) && isscalar(val.mv) && ...
                            isfinite(val.mv) && val.mv > 0 && isreal(val.mv)
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
                            isreal(val.Nafter)
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
                    
                case 'rho_Lambda_thres'
                    if isscalar(val.rho_Lambda_thres) && isreal(val.rho_Lambda_thres)
                        obj.parameters.(fn{:}) = abs(val.(fn{:}));
                    else
                        error('parameters.rho_Lambda_thres must be a real number.');
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
                            isfinite(val.n_tunnel_max) && ...
                            mod(val.n_tunnel_max,1) == 0 && isreal(val.n_tunnel_max)
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.n_tunnel_max must be a real, finite integer.');
                    end
                    
                case 'n_recycle'
                    if isnumeric(val.n_recycle) && isscalar(val.n_recycle) && ...
                            isfinite(val.n_recycle) && val.n_recycle > 0 && ...
                            mod(val.n_recycle,1) == 0 && isreal(val.n_recycle)
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.n_recycle must be a real, finite, positive integer.');
                    end
                    
                case 'seed'
                    if isnumeric(val.seed) && isscalar(val.seed) && ...
                            isfinite(val.seed) && val.seed > 0 && ...
                            mod(val.seed,1) == 0 && isreal(val.seed)
                        obj.parameters.(fn{:}) = val.(fn{:});
                    elseif isscalar(val.seed) && isnan(val.seed)
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.seed must be a real, finite, positive integer.');
                    end
                    
                case 'outfile'
                    if ischar(val.outfile)
                        obj.parameters.outfile = val.outfile;
                    else
                        error('Output file name must be a string.');
                    end
            end
        end
    end
    
end

end

function varargout = build_potential(f,Mv,Mh,d_range)
    
    nd = min(length(f),nargout)-1;
    
    % prefactor = arrayfun(@(n) p.mv.^4 * p.mh^(-n),0:3,'Un',0);
    % potential = cellfun(@(grf,A) @(phi) A*grf(phi/mh),f,prefactor,'Un',0);
    % [V,Vp,Vpp] = deal(potential{1:3});
    
    Mv4 = Mv.^4;
    
    if nargin < 4
        for d = 0:nd
            grf = f{1+d};
            mhmd = Mh^(-d);
            varargout{d+1} = @(phi) Mv4 * mhmd * grf(phi/Mh); % reshape(grf(phi(:).'),size(phi));
        end
        for d = nd+1:max(nd,nargout-1)
            varargout{d+1} = function_handle.empty();
        end
    else
        for ii = 1:length(d_range)
            d = d_range(ii);
            grf = f{1+d};
            mhmd = Mh^(-d);
            varargout{ii} = @(phi) Mv4 * mhmd * grf(phi/Mh); % reshape(grf(phi(:).'),size(phi));
        end
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


% %             if shift_potential == 1
% %                 % Rescale mv so that Q matches observation
% %                 Q_target = 2e-5;
% %                 phiexit = phi(3);
% %                 Q = sqrt(V(phiexit)/(150*(Vp(phiexit)./V(phiexit)).^2/2)*kappa^2*obj.m_Pl^2)/pi;
% %                 mv = obj.parameters.mv*sqrt(Q_target/Q);
% %                 [V0,Vp,Vpp] = build_potential(f,mv*obj.m_Pl,mh*obj.m_Pl);
% %                 rho_Lambda = rho_Lambda*(Q_target/Q)^2;
% %                 V = @(x) V0(x) - rho_Lambda;
% %             else
