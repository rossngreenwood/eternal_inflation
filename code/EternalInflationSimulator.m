classdef EternalInflationSimulator < handle

methods (Access = public)
    
    function main(obj)
        % Collect data from inflation simulations pertaining to the
        % manifestation of false-vacuum, stochastic, and topological
        % eternal inflation, as well as CMB observables.
        
        warning('off','MATLAB:integral:MaxIntervalCountReached');
        
        p = obj.parameters;
        
        % Seed the random number generator
        if isa(p.randstream,'RandStream')
            RandStream.setGlobalStream(p.randstream);
        elseif ~isnan(p.seed)
            rng(p.seed);
        end
        
        fid = fopen(p.outfile,'wt'); % Create output text file, write metadata header
        fprintf(fid,'%E,%d,%.4G,%.4G,%.4G,%d,%.2f,%s,%d,%d,%.4G,%.3G,%d,%d,%d\r\n',...
            p.n_iter,p.cores,p.mv,p.mh,obj.m_Pl,p.kmax,p.gamma,p.measure,p.n_tunnel_max,...
            p.lambdascreenmode,p.mv_offset_max,p.Nafter,p.seed,p.n_recycle,p.flag_screen_basin);
        fclose(fid);
        
        Mv = p.mv*obj.m_Pl; % Mass scales in natural units
        Mh = p.mh*obj.m_Pl;
        
        phi0_spacing = 8*Mh; % Separation between starting points on the same potential
        phi0 = phi0_spacing*(-floor(p.n_recycle/2):floor(p.n_recycle/2)-1+mod(p.n_recycle,2));
        
        i_iter = 0;
        i_success = 0;
        while true, i_iter = i_iter + 1;
            
            if p.n_iter > 0 && i_iter > p.n_iter
                break
            elseif p.n_iter < 0 && (i_success > abs(p.n_iter) || i_iter > 1e9)
                break
            end
            if mod(i_iter,1e5) == 0
                disp(num2str(i_iter));
            end
            
            data_out = nan(1,19); % Initialize output array
            
            for goto = NaN % on break, goto record data
                
                record_flag = 0;
                
                %% Draw random potential and set initial field value
                
                ii = mod(i_iter-1,p.n_recycle)+1;
                phiinit = phi0(ii); % Initial field value
                
                % Generate dimensionless potential function f and
                % compute potential values at starting points
                if mod(i_iter,p.n_recycle) == 1
                    [ak,f{1:3}] = obj.gaussian_random_field_1D(p.kmax,p.gamma);
                    V0   = Mv^4*f{1}(phi0/Mh);
                end
                
                % Maximum amount by which to shift potential
                V_offset_max = (p.mv_offset_max*obj.m_Pl)^4;
                
                switch p.measure
                    
                    case {'A','D','E'}
                        
                        % Check if potential energy is already too negative
                        if V0(ii) < -V_offset_max
                            data_out(2) = 1; break
                        else
                            data_out(2) = 0;
                        end
                        
                    case {'B','C'}
                        
                        if mod(i_iter,p.n_recycle) == 1
                            Vp0  = Mv^4*f{2}(phi0/Mh)/Mh;
                            Vpp0 = Mv^4*f{3}(phi0/Mh)/Mh^2;
                        end
                        
                        % Check if slow roll is valid at starting point
                        % Saturate the bound to give most leeway
                        if V0(ii) < -V_offset_max
                            data_out(2) = 1; break
                        else
                            data_out(2) = 0;
                        end
                        
                        % Abort if slow roll is invalid for Measure B
                        if strcmpi(p.measure,'B') && data_out(2) > 0
                            break
                        end
                        
                end
                
                %% Screen models based on vacuum energy
                
                % Find the first local minimum
                % Shift the potential up by mv_offset_max so that
                % find_phistop doesn't give up if it goes slightly negative
                [xstop] = find_phistop(phiinit/Mh,ak,...
                    1,-V_offset_max/Mv^4,1,p.lambdascreenmode);
                
                if isnan(xstop)
                    data_out(2) = 5; break % Vacuum energy is too negative
                elseif V_offset_max > 0
                    % Check for minima with sufficiently small energy
                    V_offset = obj.find_rho_offset(ak,f{1},xstop);
                    if isnan(V_offset), break, end
                    data_out(4) = V_offset/obj.m_Pl^4;
                elseif p.lambdascreenmode && f{1}(xstop) < 0
                    data_out(2) = 5; break % Vacuum energy is negative
                else
                    V_offset = 0;
                end
                
                %% Set phistart
                
                Vstart  = []; Vpstart = [];
                
                switch p.measure
                    
                    case {'A','D','E'}
                        
                        % Look uphill for the (dim-less) potential peak
                        xpeak = find_phipeak(phiinit/Mh,ak,1,0,1);
                        
                        % Move phiinit to peak; tag with fall direction
                        phiinit = xpeak*Mh + 1i*sign(phiinit-xpeak*Mh);
                        
                        % Check for slow roll at peak (assume \epsilon_V << 1)
                        if abs(f{3}(xpeak)/(f{1}(xpeak)-V_offset/Mv^4)/Mh^2)/(8*pi/obj.m_Pl^2) > 1
                            % Require slow roll at the peak for Measure A
                            data_out(2) = 3;
                            if strcmpi(p.measure,'D')
                                try
                                    phistart = obj.find_phistart_attractor(xpeak*Mh,@(x) Mv^4*f{1}(x/Mh)-V_offset,...
                                        @(x) Mv^4*f{2}(x/Mh)/Mh,@(x) Mv^4*f{3}(x/Mh)/Mh^2,Mh);
                                catch me
                                    disp(me.message);
                                    break
                                end
                            elseif strcmpi(p.measure,'E')
                                 [phistart,~] = obj.find_phistart(xpeak*Mh,@(x) Mv^4*f{1}(x/Mh)-V_offset,@(x) Mv^4*f{2}(x/Mh)/Mh,@(x) Mv^4*f{3}(x/Mh)/Mh^2,Mh,obj.m_Pl);
                                 phiinit = phistart;
                            else
                                break
                            end
                        else
                            phistart = xpeak*Mh; % Start at maximum
                        end
                        
                    case {'B','C'}
                        
                        if data_out(2) == 0 && V_offset == 0
                            phistart = phiinit; % phiinit already inflates
                            Vstart   = V0(ii)-V_offset;
                            Vpstart  = Vp0(ii);
                        elseif (f{2}(phiinit/Mh)/(f{1}(phiinit/Mh)-V_offset/Mv^4)/Mh)^2/(16*pi/obj.m_Pl^2) < 1 && ...
                                abs(f{3}(phiinit/Mh)/(f{1}(phiinit/Mh)-V_offset/Mv^4)/Mh^2)/(8*pi/obj.m_Pl^2) < 1
                            phistart = phiinit; % phiinit inflates after shift
                        elseif strcmpi(p.measure,'C')
                            % Look for any inflation downhill from phiinit
                            % phistart = obj.find_phistart(phiinit/Mh,@(x) f{1}(x)-V_offset/Mv^4,...
                            %     @(x) f{2}(x)/Mh,@(x) f{3}(x)/Mh^2,1,obj.m_Pl)*Mh;
                            phistart = obj.find_phistart_attractor(phiinit,@(x) Mv^4*f{1}(x/Mh)-V_offset,...
                                @(x) Mv^4*f{2}(x/Mh)/Mh,@(x) Mv^4*f{3}(x/Mh)/Mh^2,Mh);
                        else
                            break
                        end
                        
                end
                
                %% Simulate inflation with tunneling events
                
                record_flag = 1;
                
                phi    = nan(1+p.n_tunnel_max,6);
                status = nan(1+p.n_tunnel_max,1);
                Ntotal = nan(1+p.n_tunnel_max,1);
                
                last_valid = nan;
                weight_E = 1;
                
                for i_tunnel = 0:p.n_tunnel_max, it = 1 + i_tunnel;
                    
                    if i_tunnel > 0
                        
                        %% Check for false-vacuum tunneling
                        
                        % Look for an instanton tunneling solution
                        [phitunnel,log_tunnel_rate,flag_hawking_moss,flag_fv_eternal] = ...
                            obj.check_false_vacuum_eternal(ak,phi(i_tunnel,6),V,Vp,Vpp,...
                            V_offset,p.n_tunnel_max+1-i_tunnel);
                        
                        if isnan(phitunnel), break, end
                        if V(phi(i_tunnel,6)) <= 0, break, end
                        
                        %% Handle false-vacuum tunneling
                        
                        record_flag = 2;
                        
%                         data_out(5)  = max(0,data_out(5))  || flag_fv_eternal;   % Any false-vacuum eternal inflation?
                        data_out(18) = max(0,data_out(18)) || flag_hawking_moss; % Any Hawking-Moss instanton?
                        if ~(data_out(6) <= log_tunnel_rate)
                            data_out(19) = V(phi(i_tunnel,6))/obj.m_Pl^4;        % Energy of the false vacuum
                        end
                        data_out(6)  = min(data_out(6),log_tunnel_rate);         % Smallest tunneling rate
                        
                        if flag_hawking_moss
                            phiinit = phitunnel + 1i*sign(phitunnel-phi(i_tunnel,6));
                        else
                            phiinit = phitunnel;
                        end
                        
                        if (Vp(phitunnel)/V(phitunnel))^2/(16*pi/obj.m_Pl^2) < 1 && ...
                                abs(Vpp(phitunnel)/V(phitunnel))/(8*pi/obj.m_Pl^2) < 1
                            phistart = phitunnel; % phitunnel inflates
                        else
                            % Look for any inflation downhill from phiinit
                            phistart = obj.find_phistart(phitunnel,V,Vp,Vpp,Mh,obj.m_Pl);
                        end
                        
                        Vstart  = []; Vpstart = [];
                        
                    end
                    
                    %% Detect slowroll
                    
                    if isnan(phistart)
                        % SRA never found to be valid, even near minimum
                        status(it,1) = 4;
                        Ntotal(it,1) = NaN;
                        if i_tunnel == 0
                            [V0,Vp,Vpp] = obj.build_potential(f,Mv,Mh);
                            V = @(x) V0(x) - V_offset;
                        end
                    else
                        [phi(it,:),status(it,1),Ntotal(it,1),V,Vp,Vpp] = ...
                            obj.simulate_slowroll(f,phistart,V_offset,...
                            ~isreal(phiinit) && phistart == real(phiinit),Vstart,Vpstart);
                    end
                    
                    if ~isnan(phi(4)) && sign(phi(4)-phi(3)) ~= sign(phi(5)-phi(4))
                        break
                    end
                    
                    phi(it,2) = phiinit;
                    
                    % Set phi(:,6) = phistop, location of the local minimum
                    if status(it,1) == 4 && isfinite(phistart)
                        % Already set in simulate_slowroll
                    elseif i_tunnel == 0
                        phi(it,6) = xstop*Mh; % Found above
                    elseif isfinite(phistart)
                        phi(it,6) = find_phistop(phistart/Mh,...
                            ak,1,V_offset/Mv^4,1,0)*Mh;
                    else
                        phi(it,6) = find_phistop(real(phiinit)/Mh,...
                            ak,1,V_offset/Mv^4,1,0)*Mh;
                    end
                    
                    % For Measure E only:
                    if strcmpi(p.measure,'E') && i_tunnel == 0 && ~isnan(phi(1,4))
                        phiexit = phi(1,4);
                        % Report ratio of size of interval where N>55
                        % is possible to size of half-basin
                        weight_E = abs((real(phiinit)-phiexit)/(real(phiinit)-phi(1,6)));
                        % Resample phistart between phiexit and old
                        % phistart
                        if phistart ~= real(phiinit) % Inflation starts below maximum
                            tmp = phiexit + (real(phiinit)-phiexit)*rand();
                            if sign(phistart-phiexit) == sign(phistart-tmp)
                                % Only change phistart if the new value
                                % is lower on the slope than the old
                                phistart = tmp;
                            end
                        else % Inflation starts at maximum
                            phistart = phiexit + (phistart-phiexit)*rand();
                        end
                        phi(1,3) = phistart;
                    end
                    
                    if i_tunnel > 0 && (isnan(phi(it,6)) || V(phi(it,6)) < 0)
                        break % Vacuum energy is negative
                    elseif ~isnan(phi(it,5))
                        if sign(phi(it,5)-phi(it,3)) == sign(phi(it,6)-phi(it,5))
                            last_valid = it; % Inflation ends - OK
                        else
                            % find_phiend didn't detect local min
                            status(it,1) = 4;
                            Ntotal(it,1) = NaN;
                            phi(it,[4 5]) = [NaN, NaN];
                        end
                    end
                    
                end
                
                if ~isnan(phi(1,4)) && sign(phi(1,4)-phi(1,3)) ~= sign(phi(1,5)-phi(1,4))
                    break
                end
                    
                %% Screen models based on amount of inflation
                
                % Only fv eternal if inflation ends in another basin
                data_out(5) = (last_valid > 1);
                
                if isnan(last_valid)
                    % Inflation doesn't end in any basin?
                    last_valid = find(isfinite(status),1,'last');
                    if isempty(last_valid), break, end
                end
                
                data_out(2) = status(last_valid);
                data_out(3) = Ntotal(last_valid);
                
                phi    = phi(1:last_valid,:);
                Ntotal = Ntotal(1:last_valid);
                
                if status(last_valid,1) == 4
                    break % No end of inflation
                elseif ~(Ntotal(last_valid) > p.Nafter)
                    break % Not enough e-foldings
                end
                
                %% Fill in gaps
                
                % Compute e-folds precisely if haven't already
                if isfinite(Ntotal(end)) && Ntotal(end) > 1.3*obj.parameters.Nafter
                    kappa = 8*pi/obj.m_Pl^2;
                    dlna_dphi = @(phi) (-kappa*V(phi)./Vp(phi));
                    Ntotal(end) = integral(dlna_dphi,phi(end,3),phi(end,5));
                    data_out(3) = Ntotal(end);
                    if ~(Ntotal(end) > p.Nafter)
                        break
                    end
                end
                
                record_flag = 3;
                if p.n_iter < 0
                    i_success = i_success + 1;
                end
                
                % Set phi(:,1) = phipeak
                phi(:,1) = real(phi(:,2));
                for it = 1:size(phi,1)
                    if ~isreal(phi(it,2)) && it == 1, continue, end
                    phi(it,1) = find_phipeak(real(phi(it,2))/Mh,...
                        ak,1,V_offset/Mv^4,1)*Mh;
                end
                
                %% Check for stochastic and topological eternal inflation
                
                for b = 0:size(phi,1)-1 % Loop over basins
                    
                    % Stochastic eternal inflation
                    [numStochEpochs,NSinceStoch,NStochastic,contig_with_max,fractal_dim] = ...
                        obj.check_stochastic_eternal(V,Vp,Vpp,phi(1+b,:));
                    data_out(7) = max(0,data_out(7)) + max(0,numStochEpochs);
                    data_out(8) = NSinceStoch; % Only keeps value from last basin
                    if ~isnan(NStochastic)
	                data_out(8) = data_out(8) + imag(NStochastic)*1i;
                        Ntotal(1+b) = real(NStochastic) + NSinceStoch + p.Nafter;
                    end
                    data_out(20) = fractal_dim;
                    
                    % Topological eternal inflation
                    flag_topological_eternal = obj.check_topological_eternal(V,Vp,Vpp,phi(1+b,2),phi(1+b,1));
                    if strcmpi(p.measure,'B') && ~isnan(contig_with_max) && contig_with_max
                        flag_topological_eternal = true;
                    end
                    data_out(9) = max(0,data_out(9)) + max(0,flag_topological_eternal);
                    
                end
                
                data_out(3) = Ntotal(end);
                
                %% Compute observables
                
                f{4} = function_handle.empty; % Construct V'''(phi)
                [~,f{1:4}] = obj.gaussian_random_field_1D(p.kmax,p.gamma,ak);
                Vppp = obj.build_potential(f,p.mv*obj.m_Pl,p.mh*obj.m_Pl,3);
                
                Nbefore = Ntotal(end) - p.Nafter; % e-folds before horizon crossing
                observables = obj.compute_observables(V,Vp,Vpp,Vppp,phi(end,:),Nbefore,obj.m_Pl);
                
                data_out(10) = observables.Q;
                data_out(11) = observables.r;
                data_out(12) = observables.n_s;
                data_out(13) = observables.alpha;
                data_out(14) = observables.n_t;
                data_out(15) = observables.dlgrho;
                data_out(16) = observables.lgOk;
                data_out(17) = observables.rho_Lambda/obj.m_Pl^4;
%                 data_out(20) = weight_E;
                
            end % for goto
            
            %% Record results
            
            obj.write_output(record_flag,p.outfile,data_out);
            
        end % while
        
    end
    
end

methods (Access = protected)
    
    function [phi,status,Ntotal,V,Vp,Vpp] = simulate_slowroll(obj,f,phistart,rho_offset,flag_start_max,Vinit,Vpinit)
        
        phi    = nan(1,6); % [phipeak, phiinit, phistart, phiexit, phiend, phistop]
        
        mv = obj.parameters.mv;
        mh = obj.parameters.mh;
        
        Ntotal = nan;
        
        kappa = 8*pi/obj.m_Pl^2;
        
        if nargin < 4 || isempty(rho_offset) || rho_offset == 0
            [V,Vp,Vpp] = obj.build_potential(f,mv*obj.m_Pl,mh*obj.m_Pl);
        else
            [V0,Vp,Vpp] = obj.build_potential(f,mv*obj.m_Pl,mh*obj.m_Pl);
            V = @(x) V0(x) - rho_offset;
        end
        
        if nargin < 5 || isempty(flag_start_max), flag_start_max = 0; end
        
        if nargin < 6 || isempty(Vinit),  Vinit  = V(phistart);  end
        if nargin < 7 || isempty(Vpinit), Vpinit = Vp(phistart); end
        
        %% Identify phistart
        
        phi(3) = phistart;
            
        %% Identify phiend, phiexit, and phistop
        
        % Find the end of slow roll inflation
        [phiend,status] = obj.find_phiend(phistart,V,Vp,Vpp,...
            mh*obj.m_Pl,obj.m_Pl,Vinit,Vpinit);
        if status == 4 % Slow roll valid until minimum; no phiend
            phi(6) = phiend;
            return
        else
            phi(5) = phiend;
        end
        
        % Compute total number of slow roll e-folds
        dlna_dphi = @(phi) (-kappa*V(phi)./Vp(phi));
        if flag_start_max
            Ntotal = Inf; % Started at a maximum
        else
            points = linspace(phistart,phiend,max(10,abs(phistart-phiend)/mh/obj.m_Pl*10));
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
            phi(4) = phiexit(3);
            
        end
        
    end
    
    function [phitunnel,log_tunnel_rate,flag_hm,flag_eternal] = check_false_vacuum_eternal(...
            obj,ak,phifv,V,Vp,Vpp,rho_offset,n_tunnel,Vfalse)
        
        if nargin < 7 || isempty(rho_offset),  rho_offset = 0;      end
        if nargin < 8 || isempty(n_tunnel),    n_tunnel = 1;        end
        if nargin < 9 || isempty(Vfalse),      Vfalse = V(phifv);   end
        
        hbar = 1;
        log_stable_rate_cutoff = -inf; %-1e7;
        
        phitunnel = nan(1,2);
        flag_hm   = false(1,2);
        B         = inf(1,2);
        
        kappa = 8*pi/obj.m_Pl^2;
        
        %% Check for tunneling to nearby minima
        
        for lr = 1:2 % Left and right neighbor basins
            
            % Check if inflation can be successful in neighboring basins
            [any_viable,phitv,phipeak] = obj.screen_basin(...
                lr,V,Vp,Vpp,ak,rho_offset,n_tunnel,phifv,Vfalse);
            if isnan(phipeak), continue, end
%             if abs(V(phitv)) > (obj.parameters.mv*obj.m_Pl)^4/100
%                 continue
%             end
            
            mh_eff = 1/sqrt(abs(Vpp(phipeak(1))/V(phipeak(1))))/obj.m_Pl;
            
            if obj.parameters.flag_screen_basin && mh_eff > 1
                % There's almost never a CdL solution for mh > 1
                
                %% Pre-compute Hawking-Moss tunneling rate
                %  Check if a Hawking-Moss transition can
                %  give us inflation near the local maximum
                
                w_top = abs(kappa/3*V(phipeak(1)))^0.5;
                R = pi/w_top/2; Y = [phipeak(1),0,1/w_top,0,-w_top];
                B_HM = FalseVacuumInstanton.find_tunneling_suppression_static(...
                    3,kappa,V,phifv(1),R,Y);
                
                flag_hm(lr) = true;  % Hawking-Moss instanton?
                phitunnel(lr) = phipeak(1); % Field value at center of bubble
                % Tunneling suppression B = -log(\lambda) = S_bubble - S_bkgd
                B(lr) = B_HM;
                
            elseif mh_eff > 0.1 || any_viable
                
                try
                    
                    % Initialize instanton solver
                    fvi = FalseVacuumInstanton(...
                        'V',            V,...
                        'dV',           Vp,...
                        'd2V',          Vpp,...
                        'm_Pl',         obj.m_Pl,...
                        'phi_metaMin',  phifv,...
                        'phi_absMin',   phitv,...
                        'phi_bar_top',  phipeak);
                    
                    % Solve for instanton profile
                    [R,profile] = fvi.find_profile([],...
                        1e-4,... % xtol
                        1e-4,... % phitol
                        1e-2 );  % thinCutoff
                    
                catch me
                    switch me.identifier
                        case 'FalseVacuumInstanton:StableFalseVacuum'
                            continue % No tunneling
                        case 'FalseVacuumInstanton:IntegralDiverged'
                            continue % Integration failed; assume no tunneling
                        otherwise
%                             % Write potential coefs, phi_meatMin, phi_absMin to log
%                             % file, later for debugging
%                             disp(me.message);
%                             fid = fopen(obj.logfile,'at');
%                             fprintf(fid,'%s\r\n',me.message);
%                             fprintf(fid,'n_tunnel_remaining: %d\r\n',n_tunnel);
%                             fprintf(fid,'rho_offset: %.4G\r\n',rho_offset);
%                             fprintf(fid,'phi_metaMin, phi_absMin: %.4G,%.4G; ak:\r\n',[phifv,phitv]);
%                             for ii = 1:size(ak,1)
%                                 fprintf(fid,'%.4f %.4f;\r\n',ak(ii,:));
%                             end
%                             fclose(fid);
                            continue
                    end
                end
                
%                 if ~isscalar(R)
%                     x = linspace(phitv,phifv,1000); plot(x,V(x)); hold on; plot(profile(1),V(profile(1)),'o'); plot(profile(end,1),V(profile(end,1)),'x');
%                     hold off
%                     pause(2)
%                 end
                if isfinite(profile(1,1))
                    flag_hm(lr) = isscalar(R);  % Hawking-Moss instanton?
                    phitunnel(lr) = profile(1,1); % Field value at center of bubble
                    % Tunneling suppression B = -log(\lambda) = S_bubble - S_bkgd
                    B(lr) = fvi.find_tunneling_suppression(R,profile);
                end
                
            else
                continue
            end
            
        end
        
        % Handle no-tunneling case
        if all(isnan(phitunnel))
            phitunnel       = NaN;
            log_tunnel_rate = NaN;
            flag_hm         = false;
            flag_eternal    = false;
            return
        end
        
        %% Choose bubble with the larger transition rate
        
        [log_tunnel_rate,imax] = max(-B/hbar); % B has units of kg*m
        
        % Find the new value of phi after tunneling
        if log_tunnel_rate > log_stable_rate_cutoff
            phitunnel = phitunnel(imax);
            flag_hm = flag_hm(imax);
        else
            phitunnel = NaN;
            flag_hm = false;
        end
        
        % Determine if inflation is eternal; H4 = (kappa/3*V(phistop))^2;
        flag_eternal = ( log_tunnel_rate < log(9/4/pi) );
        
    end
    
    function [numStochEpochs,NSinceStoch,NStochastic,contig_with_max,fractal_dim] = check_stochastic_eternal(obj,V,Vp,Vpp,phi)
        
        p = obj.parameters;
        
        kappa = 8*pi/obj.m_Pl^2;
        hbar  = 1;
        
        numStochEpochs = nan; % Number of intervals of eternal inflation
        NSinceStoch    = nan; % Number of e-foldings between SEI breakdown and exit scale
        NStochastic    = nan; % Expected number of e-folds 
        contig_with_max = nan;
        fractal_dim    = nan;
        
        phiinit  = phi(2);
        phistart = phi(3);
        phiexit  = phi(4);
        phiend   = phi(5);
        
        % Initialized at peak and the peak inflates
        flag_starts_at_peak = (imag(phi(2)) ~= 0) && (phistart == real(phi(2)));
        
        phiscale = p.mh*obj.m_Pl;
        
        %% Find intervals of stochastic eternal inflation
        
        Vstart  = V(phistart);
        Vpstart = Vp(phistart);
        
        sgn    = sign(Vp(phistart));
        dphi   = -0.001*sgn*max(0.01*phiscale,min(phiscale^2 * abs(Vpstart./Vstart),phiscale));
        
        phibreak = nan; % Field values where SEIC changes sign
        off2on   = nan; % True if SEIC goes from (+) to (-) moving downhill
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
        
        if isnan(phibreak(1))
            % Reached phiend without satisfying 
            % eternal inflation criterion
            return
        end
        
        %% Find field values at breakdown precisely using binary search
        
        dphi = 1e-7;
        
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
        
        %% Screen phistart and phiend
        
        if off2on(1) == 0 % SEI is satisfied at phistart
            % Add phistart as the past boundary of the first interval
            phibreak = [phistart phibreak];
            off2on   = [1        off2on];
            % Determine if stochastic interval is contiguous with maximum
            phipeak = phi(1);
            dphi   = 0.001*sgn*max(0.01*phiscale,min(phiscale^2 * abs(Vpstart./Vstart),phiscale));
            ind = 1; ii = 1;
            phisearch = phistart*ones(1,batch);
            contig_with_max = true;
            while sgn*(phisearch(ii)-phipeak) < 0
                if mod(ind,batch) == 1
                    phi_last = phisearch(end);
                    phisearch = phi_last + cumsum(dphi*step.^(ind-1:ind+batch-2));
                    Vend   = V(phisearch);
                    Vpend  = Vp(phisearch);
                end
                ii = mod(ind-1,batch)+1;
                if (hbar^0.5*(kappa*Vend(ii))^1.5 - 2*pi*sqrt(3)*(0.607)*abs(Vpend(ii))) < 0
                    contig_with_max = false;
                    break
                end
                ind = ind + 1;
            end
        end
        
        if off2on(end) == 1 % SEI is satisfied at phiend
            % Require that SEI ends when inflation ends, if not sooner
            phibreak(end+1) = phiend;
            off2on(end+1)   = 0;
        end
        
        %% Compute NStochastic
        
        NStochastic = 0;
        ii = 1;
        while true
            if ii > length(phibreak)    
                break
            end
            DeltaPhi = abs(phibreak(ii+1)-phibreak(ii));
            deltaPhi = sqrt(kappa*V((phibreak(ii+1)+phibreak(ii))/2)/3)/2/pi;
            if phibreak(ii) == phistart && (phistart == real(phiinit)) && (DeltaPhi/deltaPhi) < 1
                % fluctuations are larger than fluctuation-dominated
                % interval near the maximum; don't count it
                phibreak(ii:ii+1) = [];
                off2on(ii:ii+1)   = [];
                NStochastic = real(NStochastic) + (DeltaPhi/deltaPhi)*1i
                continue
            end
            NStochastic = NStochastic + (DeltaPhi/deltaPhi)^2;
            ii = ii+2;
        end
        
        if isempty(phibreak)
            return
        end
        
        numStochEpochs = nnz(~off2on);
        
        %% Compute fractal dimension
        
        if phibreak(1) == phistart && (phistart == real(phiinit))
            H = sqrt(kappa*V(phibreak(1))/3);
            M = -Vpp(phistart);
            fractal_dim = -M^2/3/H^2;
        end
        
        %% Compute # of e-folds after past SEI breakdown and before phiexit
        
        dlna_dphi = @(phi) -V(phi)./Vp(phi); % Integrate this to find N_e
        
        if isnan(phiexit)
        elseif (hbar^0.5*kappa*V(phiexit)).^1.5 - 2*pi*sqrt(3)*(0.607)*abs(Vp(phiexit)) > 0
            NSinceStoch = 0; % Eternal at phiexit
        elseif any((phibreak(~off2on)-phiexit)*Vp(phiexit) > 0)
            % Closest breakdown point of SEI to phiexit
            phibreak_last = phibreak((phibreak-phiexit)*Vp(phiexit) > 0 & ~off2on);
            [~,imin] = min(abs(phibreak_last-phiexit));
            % Integrate e-foldings
            NSinceStoch = integral(@(x) dlna_dphi(x),phibreak_last(imin),phiexit);
        end
        
    end
    
    function [flag_topological_eternal] = check_topological_eternal(obj,V,Vp,Vpp,phistart,phipeak)
        
        p = obj.parameters;
        
        flag_topological_eternal = nan;
        
        %% Topological Eternal Inflation
        
        kappa = 8*pi/obj.m_Pl^2;
        hbar  = 1;
        
        % Determine whether quantum fluctuations could result in at least
        % one Hubble volume descending toward a different minimum of V(phi)
        if isreal(phistart) % Not starting at maximum
            H = sqrt(kappa/3*V(phistart));
            dphi = hbar^0.5*H/2/pi; % Amplitude of quantum fluctuation: P_\phi(k) = (dphi)^2 = \hbar * (H/2/\pi)^2
            Dphi = -Vp(phistart)/(3*H^2); % Classical field excursion
            if 1/2*erfc(abs(phistart+Dphi-phipeak)/(sqrt(2)*dphi)) < exp(-3)
                % Not close enough to maximum; assume no domain wall forms
                return
            end
        end
        
        %% Assuming a domain wall forms around maximum, does it persist?
        
        if abs(Vpp(phipeak)/V(phipeak))/(8*pi/obj.m_Pl^2) > 1
            % No inflation at maximum
            flag_topological_eternal = false;
        else
            flag_topological_eternal = true;
        end
        
        return
        
        %%
        
        % Find value of phi at the domain wall boundary 
        % where inflation breaks down
        [phiedge,status] = obj.find_phiend(phipeak,V,Vp,Vpp,p.mh*obj.m_Pl,obj.m_Pl,[],[],false,false);
        if status == 1
            % Inflation ends when the potential goes negative
            flag_topological_eternal = false;
            return
        elseif status == 4
            % Inflation doesn't end
            flag_topological_eternal = true;
            return
        end
        
        % Check whether the amplitude of fluctuations at the maximum places
        % more than (1-1/e) outside the inflating interval
        if ~isreal(phistart) % Not starting at maximum
            H = sqrt(kappa/3*V(real(phistart)));
            dphi = hbar^0.5*H/2/pi; % Amplitude of quantum fluctuation: P_\phi(k) = (dphi)^2 = \hbar * (H/2/\pi)^2
        end
        if erf(abs(phiedge-phipeak)/(sqrt(2)*dphi)) < exp(-1)
            % Fluctuations are too big at peak
            flag_topological_eternal = false;
            return
        end
        
        % Find value of phi that will descend to phiedge in time <= t_H
        dphi_dlna = @(lna,phi) -Vp(phi)./V(phi)/kappa;
        [~,phiout] = ode23(dphi_dlna,[1,0.5,0],phiedge);
        phistar = phiout(3);
        
        if isnan(phistar)
            return
        end
        
        % If the expansion of the inner portion of the domain wall
        % where (phipeak < phi < phistar) replaces loss of the outer
        % wall where (phistar < phi < phiedge), then inflation is eternal
        % (Problem is effectively 1D - the width of domain wall.)
        if abs(phistar-phipeak) > abs(phiedge-phipeak)*exp(-1)
            flag_topological_eternal = true;
        end
        
    end
    
    function [rho_offset,status] = find_rho_offset(obj,ak,f,xstop)
        
        p = obj.parameters;
        
        % Mass scales in natural units
        Mv = p.mv*obj.m_Pl;
        
        % Check if there is a minimum close to mv_offset_max
        
        if isnan(xstop)
            % Vacuum energy is too negative; abort
            status = 5;
            rho_offset = nan;
            return
        end
        
        f_offset_max = p.mv_offset_max^4/p.mv^4;
        
        fstop = f(xstop);
        if abs(fstop) < f_offset_max
            % Vacuum energy is already close to a possible
            % value for lambda in the starting basin
            rho_offset = fstop*Mv^4;
        else
            
            % Check if a neighboring vacuum is close to the
            % observed value for Lambda
            
            phitv   = {nan(1,p.n_tunnel_max),nan(1,p.n_tunnel_max)};
            phipeak = {nan(1,p.n_tunnel_max),nan(1,p.n_tunnel_max)}; % Not using this
            phifv   = repmat({[xstop nan(1,p.n_tunnel_max-1)]},1,2);
            
            rho_offset = nan(1,2);
            rho_basin  = nan(1,2);
            
            for lr = 1:2
                for i_tunnel = 1:p.n_tunnel_max
                    [phitv{lr}(i_tunnel),phipeak{lr}(i_tunnel)] = find_phinextmin(...
                        phifv{lr}(i_tunnel),ak,1,...
                        -f_offset_max,...
                        1,1,2*(lr-1)-1,1);
                    if ...
                            isnan(phitv{lr}(i_tunnel)) || ...
                            f(phitv{lr}(i_tunnel)) < 0 || ...
                            diff(f([phitv{lr}(i_tunnel),phifv{lr}(i_tunnel)])) < 0
                        break
                    elseif abs(f(phitv{lr}(i_tunnel))) < f_offset_max
                        rho_offset(lr) = f(phitv{lr}(i_tunnel))*Mv^4;
                        rho_basin(lr)  = i_tunnel;
                        break
                    elseif i_tunnel < p.n_tunnel_max
                        phifv{lr}(i_tunnel+1) = phitv{lr}(i_tunnel);
                    end
                end
            end
            
            if any(~isnan(rho_offset))
                % Choose the offset that zeros out the closest
                % basin or that involves the smallest shift
                if rho_basin(1) == rho_basin(2)
                    [~,ia] = sort(abs(rho_offset));
                else
                    [~,ia] = sort(rho_basin); % Always chooses non-NaN
                end
                rho_offset = rho_offset(ia(1));
            else
                rho_offset = nan;
                return
            end
            
        end
        
    end
    
    function [any_viable,next_minimum,next_peak,phifv,Vfv] = screen_basin(...
            obj,direction,V,Vp,Vpp,ak,rho_offset,n_tunnel_remain,phistop,Vstop)
        % Look ahead up to next n_tunnel_remaining adjacent basins. If we
        % can't get enough e-foldings on the other side the potential
        % barrier, don't bother computing the tunneling rate
        
        Mh = obj.parameters.mh*obj.m_Pl;
        
        kappa = 8*pi/obj.m_Pl^2;
        log_stable_rate_cutoff = -1e2;
        
        near_minima = NaN;
        
        phifv   = phistop;
        phipeak = nan;
        Vfv     = Vstop;
        
        for b = 1:n_tunnel_remain
            
            %% Find next neighboring minimum
            
            [xmin,xpeak] = find_phinextmin(phifv(b)/Mh,...
                ak,1,rho_offset/obj.parameters.mv^4/obj.m_Pl^4,1,1,2*(direction-1)-1,1);
            near_minima(b) = xmin*Mh;
            phipeak(b)     = xpeak*Mh;
            
            if isnan(near_minima(b))
                any_viable = false; break % Next minimum is negative
            elseif V(near_minima(b)) > Vfv(b)
                any_viable = false; break % "True" vacuum > "false" vacuum
            end
            
            if ~obj.parameters.flag_screen_basin
                % Just find next minimum
                any_viable = true; break
            end
            
            %% Find barrier edge location
            % Use this starting point to compute the maximum amount of
            % inflation that could occur in the next potential basin
            % Vectorized binary search
            
            phi_tol = abs(phifv(b)-near_minima(b))*1e-10;
            
            phimin = phipeak(b);
            phimax = near_minima(b);
            sgn = sign(phimax-phimin);
            
            if phimin > phimax
                phitemp = phimin;
                phimin = phimax;
                phimax = phitemp;
            end
            phibar = 0.5*(phimin + phimax);
            phisep = abs(phimax-phimin);
            
            Vtarget = 0.05*V(near_minima(b)) + 0.95*V(phipeak(b));
%             Vtarget = V(phifv);
            
            nbits = 5;
            ind = 1;
            while abs(phisep)*2^(1-ind) > abs(phi_tol)
                if mod(ind,nbits-1) == 1
                    phisep = (phimax-phimin);
                    phi_range = (phimin + phisep*(2^-nbits)):(phisep*(2^-nbits)):(phimax - phisep*(2^-nbits));
                    fun_vals = sgn*(Vtarget-V(phi_range));
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
            
            %% Find the start of inflation in the neighboring basin
            
            Vbar = V(phibar);
            if (Vp(phibar)/Vbar).^2/(2*kappa) > 1
                status_start = 1;
            elseif abs(Vpp(phibar)/Vbar)/(kappa) > 1
                status_start = 2;
            else
                status_start = 0;
            end
            
            if status_start == 0 % Inflation at barrier edge
                phistart = phibar;
            else
                % Look for inflation farther down the slope
                [phistart,status_start] = obj.find_phistart(...
                    phibar,V,Vp,Vpp,Mh,obj.m_Pl);
            end
            
            %% Compute maximum inflation starting from barrier edge
            
            if status_start == 0 % Found inflation (should always happen)
                
                % Find the end of slow roll inflation in neighbor basin
                [phiend,status_end] = obj.find_phiend(phistart,V,Vp,Vpp,...
                    Mh,obj.m_Pl,[],[],false,(b < n_tunnel_remain));
                
                if status_end ~= 4 % Inflation ends in this basin
                    % Compute total number of slow roll e-folds
                    dlna_dphi = @(phi) (-kappa*V(phi)./Vp(phi));
                    points = linspace(phistart,phiend,max(10,abs(phistart-phiend)/Mh*10));
                    Ntotal = trapz(points,dlna_dphi(points));
                    if  Ntotal > 0.7*obj.parameters.Nafter
                        Ntotal = integral(dlna_dphi,phistart,phiend);
                    end
                end
                
            end
            
            %% Check if there is enough inflation
            
            if status_start == 0 && status_end == 4
                % No end of inflation; not viable
            elseif status_start == 0 && Ntotal >= obj.parameters.Nafter
                % Tunneling basin is viable
                any_viable = true; break
            end
%             elseif abs(Vpp(phipeak(b))/V(phipeak(b)))/(kappa) < 1
%                 % Second derivative slow roll conditions met at peak
%                 
%                 %% Pre-compute Hawking-Moss tunneling rate
%                 %  Check if a Hawking-Moss transition can
%                 %  give us inflation near the local maximum
%                 
%                 w_top = abs(kappa/3*V(phipeak(b)))^0.5;
%                 R = pi/w_top/2; Y = [phipeak(b),0,1/w_top,0,-w_top];
%                 B_HM = FalseVacuumInstanton.find_tunneling_suppression_static(...
%                     3,kappa,V,phifv(b),R,Y);
%                 
%                 if -B_HM >= log_stable_rate_cutoff
%                     % Hawking-Moss instanton is viable
%                     any_viable = true; break
%                 end
%                 
%             end
            
            if b < n_tunnel_remain
                phifv(b+1) = near_minima(b);
                Vfv(b+1)   = V(phifv(b+1));
            else
                any_viable = false; % No more chances to tunnel
            end
            
        end % If loop exits naturally, no viable basin was found
        
        next_minimum = near_minima(1);
        next_peak    = phipeak(1);
        
    end
    
    function [phistart] = find_phistart_attractor(obj,phipeak,V,Vp,Vpp,phiscale)
        
        %% Find start of slow roll below the maximum
        
        [phistart,status] = obj.find_phistart(phipeak,V,Vp,Vpp,phiscale,obj.m_Pl);
        
        if isnan(phistart)
            return
        end
        
        %% Test the maximum value |dphi/dt| = sqrt(2*V(phistart))
        
        % Positive means we're heading down
        ysign = sign(phistart-phipeak);
        Vstart = V(phistart);
        
        y0 = [phistart, ysign*sqrt(2*Vstart), 1, sqrt(V(phistart)/3/obj.M_Pl^2)];
        H = sqrt(V(phipeak)/3/obj.M_Pl^2);
        y0 = [phipeak+ysign*(H/2/pi), 0, 1, H]; % Start 1-\sigma away from peak
        
        % dY is the ODE that we use
        dY = @(r,y) obj.equation_of_motion(V,Vp,r,y,obj.m_Pl).';
        
        ie = [];
        tmin = abs(y0(1)/y0(2))/10;
        count = 0;
        while isempty(ie)
            
            options = odeset(...
                'Events',   @(~,y) obj.ode_events(y,ysign,phistart) );
            [r1,y1,re,~,ie] = ode23s(dY,[0,tmin],y0.',options);
            
            tmin = tmin*2;
            count = count+1;
            if count > 10
                error('EternalInflationSimulator:NoTerminatingEvents',...
                    'No terminating events were triggered.')
            end
            
        end
        
        %% Return phistart = nan if phi picks up too much kinetic energy
        
        switch ie(end)
            case 1 % Passed the peak
                % Means kinetic energy is within attractor at phistart
                if abs(y1(2)) > sqrt(2*Vstart)
                    phistart = nan;
                end
%             case 2 % Didn't make it to the peak
%                 phistart = nan;
            otherwise
                error('FalseVacuumInstanton:IntegralDiverged',...
                    'ODE solver failed to integrate instanton solution.');
        end
        
    end
    
end

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
    
    function [varargout] = build_potential(f,Mv,Mh,d_range)
        
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
    
    function [phistart,status] = find_phistart(phiinit,V,Vp,Vpp,phiscale,m_Pl)
        % Find the value of phi at the end of slow roll inflation
        
        sgn    = sign(Vp(phiinit));
        phi1 = phiscale^2 * abs(Vp(phiinit)./V(phiinit));
        dphi   = -0.001*sgn*max(0.01*phiscale,min(phi1,phiscale));
        
        kappa = 8*pi/m_Pl^2;
        
        % Determine if inflation is valid at phiinit
        % If we start in a region with inflation, allow the field to exit
        % that region and then look for the next start of inflation
        if (Vp(phiinit)/V(phiinit)).^2/(2*kappa) > 1
            status = 2;
        elseif abs(Vpp(phiinit)/V(phiinit))/(kappa) > 1
            status = 3;
        else
            status = 0;
        end
        flag_off2on = (status > 0);
        
        step = 10^(1/16);
        ind = 1;
        batch = 40;
        phi = phiinit*ones(1,batch);
        while (~flag_off2on || status ~= 0)
            status_last = status;
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
                break
            elseif (Vpend(ii)/Vend(ii)).^2/(2*kappa) > 1
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
            case 2, fun = @(x) (Vp(x)./V(x)).^2/(2*kappa) - 1;  % eps = 1
            case 3, fun = @(x) abs(Vpp(x)./V(x))/kappa - 1;     % eta = 1
        end
        
        if ii > 1 % fun(phi1) > 0
            phi1 = phi(ii-1);
        else
            phi1 = phi_last;
        end
        phi2 = phi(ii); % fun(phi2) < 0
        
        % Make sure SRA is valid at phistart
        phistart = 0.5*(phi1 + phi2);
        while abs(phi2-phi1) > abs(dphi) || fun(phistart) > 0
            if fun(phistart) > 0
                phi1 = phistart;
            else
                phi2 = phistart;
            end
            phistart = 0.5*(phi1 + phi2);
        end
        
    end
    
    function [phiend,status] = find_phiend(phistart,V,Vp,Vpp,phiscale,m_Pl,...
            Vstart,Vpstart,lambdascreenmode,precisephistop)
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
                break
            elseif abs(Vppend(ii)/Vend(ii))/kappa > 1
                break
            end
            if sgn*Vpend(ii) < 0
                break
            end
            ind = ind + 1;
        end
        phiend = phi;
        
        if status == 0
            % Take steps until passed breakdown of slow roll
            step = 10^(1/16);                                                       %%% Tunable
            ind = 1;
            batch = 40;
            if ii == 1
                phi = phi_last*ones(1,batch);
            else
                phi = phiend(ii-1)*ones(1,batch);
            end
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
        end
        
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
        
        phiexit = phi(4);
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
        observables.dlgrho      = log10(V_exit/V(phi(5)));
        observables.lgOk        = log10(V(phi(3))/V_exit) - Nbefore*2/log(10);
        observables.rho_Lambda  = V(phi(6));
        
    end
    
    function write_output(record_flag,outfile,data_out)
        
        % Record Flags
        % 1     Potential is valid; simulated slow roll
        % 2     Checked for false vacuum tunneling
        % 3     Found a valid basin with sufficient N_e
        
        % Idx   Output          Type        Record Flag
        % 1     mv              (float,4)   1,2,3
        % 2  	status          (int)       1,2,3
        % 3     Ntotal          (float,2)   1,2,3
        % 4		rho_offset      (float,4)   3
        % 5  	flag_fv_eternal (bool)      2,3
        % 6  	log_tunnel_rate (float,4)   2,3
        % 7  	numStochEpochs  (int)       3
        % 8  	NSinceStoch     (float,2)   3
        % 9  	numTopolEpochs  (int)       3
        % 10  	Q               (float,4)   3
        % 11 	r               (float,4)   3
        % 12 	n_s             (float,4)   3
        % 13 	alpha           (float,4)   3
        % 14 	n_t             (float,4)   3
        % 15 	dlgrho          (float,4)   3
        % 16 	lgOk            (float,4)   3
        % 17 	rho_Lambda      (float,4)   3
        
        % Open output file for recording results
        % Record various amounts of data based on how far code ran
        
        switch record_flag
            case 1
                fid = fopen(outfile,'at');
                fprintf(fid,'%d,%d,%.3G,%.4G\r\n',...
                    [1 data_out(2:4)]);
                fclose(fid);
            case 2
                fid = fopen(outfile,'at');
                fprintf(fid,'%d,%d,%.3G,%.4G,%d,%.4G,%d,%.4G\r\n',...
                    [2 data_out([2:6 18 19])]);
                fclose(fid);
            case 3
                fid = fopen(outfile,'at');
                fprintf(fid,'%d,%d,%.3G,%.4G,%d,%.4G,%d,%.4G,%d,%.2G,%d,%.4G,%.4G,%.4G,%.4G,%.4G,%.4G,%.4G,%.4G,%.4G\r\n',...
                    [3 data_out([2:6 18 19 7:17 20])]);
                fclose(fid);
        end
        
    end
    
    function dy = equation_of_motion(V,dV,~,y,m_Pl)
        % Used to integrate the bubble wall.
        %
        % .. math::
        %   \frac{d^2\phi}{dr^2} 
        %     + \frac{\alpha}{\rho}\frac{d\rho}{dr}\frac{d\phi}{dr} =
        %     \frac{dV}{d\phi} \\
        %   \frac{d^2\rho}{dr^2} = 
        %     -\frac{8\pi}{3M_{Pl}^2}\left[\left(\frac{d\phi}{dr}\right)^2
        %     + V(\phi)\right]
        % 
        % Inputs
        %   y = [phi, dphi, rho, drho]
        %   r (unused)
        % Outputs
        %   dy = [dphi, d2phi, drho, d2rho]
        
        y = permute(y,circshift([1 2],[0,find(size(y) == 4)]));
        
        % Unpack variables
        phi  = y(:,1);
        dphi = y(:,2);
        a    = y(:,3);
        da   = y(:,4);
        
        kappa = 8*pi/m_Pl^2;
        
        d2phi = -dV(phi.').' - 3*da.*dphi./a;
        d2a = -kappa/3*a*(dphi.^2 - V(phi.').');
        
        dy = [dphi,d2phi,da,d2a];
        
    end
    
    function [value,isterminal,direction] = ode_events(y,ysign,phi0)
        
        phi = y(1); dphi = y(2);
        
        events = [...
            (phi-phi0)*ysign,  true,   0;...
%             dphi*ysign,        true,   0 ...
            ];
        
        value       = events(:,1);
        isterminal  = events(:,2);
        direction   = events(:,3);
        
    end
    
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
                            isfinite(val.n_iter) && val.n_iter ~= 0 && ...
                            mod(val.n_iter,1) == 0 && isreal(val.n_iter)
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.n_iter must be a positive integer.');
                    end
                    
                case 'cores'
                    if isnumeric(val.cores) && isscalar(val.cores) && ...
                            isfinite(val.cores) && val.cores ~= 0 && ...
                            mod(val.cores,1) == 0 && isreal(val.cores)
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.cores must be a positive integer.');
                    end
                    
                case 'mv'
                    
                    if isnumeric(val.mv) && isscalar(val.mv) && ...
                            isfinite(val.mv) && val.mv > 0 && isreal(val.mv)
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.mv must be a positive number.');
                    end
                    
                case 'mh'
                    
                    if isnumeric(val.mh) && isscalar(val.mh) && ...
                            isfinite(val.mh) && val.mh > 0 && isreal(val.mh)
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.mh must be a positive number.');
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
                    
                case 'mv_offset_max'
                    if isscalar(val.mv_offset_max) && isreal(val.mv_offset_max)
                        obj.parameters.(fn{:}) = abs(val.(fn{:}));
                    else
                        error('parameters.mv_offset_max must be a real number.');
                    end
                    
                case 'measure'
                    if ismember(upper(val.measure),{'A','B','C','D','E'})
                        obj.parameters.(fn{:}) = val.(fn{:});
                    else
                        error('parameters.measure must be ''A'' or ''B''');
                    end
                    
                case 'n_tunnel_max'
                    if isnumeric(val.n_tunnel_max) && isscalar(val.n_tunnel_max) && ...
                            isfinite(val.n_tunnel_max) && ...
                            mod(val.n_tunnel_max,1) == 0 && isreal(val.n_tunnel_max)
                        if val.n_tunnel_max >= 0
                            obj.parameters.(fn{:}) = val.n_tunnel_max;
                            obj.parameters.flag_screen_basin = true;
                        else
                            obj.parameters.(fn{:}) = -1*val.n_tunnel_max;
                            obj.parameters.flag_screen_basin = false;
                        end
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
                        error('parameters.seed must be a real, finite, positive integer, or NaN.');
                    end
                    
                case 'outfile'
                    if ischar(val.outfile)
                        obj.parameters.outfile = val.outfile;
                        obj.logfile = [val.outfile(1:end-4) '_log.txt'];
                    else
                        error('Output file name must be a string.');
                    end
                    
                case 'randstream'
                    if isempty(val.randstream) || isa(val.randstream,'RandStream')
                        obj.parameters.randstream = val.randstream;
                        if ~isempty(val.randstream)
                            obj.parameters.seed = val.randstream.Seed;
                        end
                    else
                        error('randstream must be an object of class RandStream.');
                    end
                    
                case 'flag_screen_basin'
                    if isscalar(val.flag_screen_basin)
                        obj.parameters.flag_screen_basin = val.flag_screen_basin;
                    else
                        error('flag_screen_basin must be a scalar logical');
                    end
            end
        end
    end
    
end

properties (SetAccess = protected)
    
    parameters = struct(...
        'seed',               NaN,...   % Seed for random number generator (integer)
        'n_iter',             1e2,...   % Number of iterations
        'mv',                 1e0,...   % Initial mass scale of the potential (multiple of m_Pl)
        'mh',                 1,...     % Mass scale of the inflaton field (multiple of m_Pl)
        'kmax',               30,...    % Largest wavenumber for GRF
        'gamma',              0,...     % Frequency dependence of GRF
        'Nafter',             55,...    % Number of e-folds between phiexit and phiend
        'lambdascreenmode',   true,...  % Throw out cases where rho_Lambda < 0?
        'measure',            'B',...   % Measure on initial conditions
        'n_tunnel_max',       1,...     % Max number of tunneling events to simulate
        'n_recycle',          4,...     % # of times to reuse same V(phi) with different phi0 value
        'mv_offset_max',      0,...     % Threshold below which vacuum energy is considered "small"
        'outfile',            '',...    % Path to output text file
        'randstream',         [],...    % Object defining random number generator
        'cores',              1,...     % Number of parallel cores on which this is running
        'flag_screen_basin',  false ...
        );
    
    logfile
    
end

properties (Constant)
    
    % Always hbar = c = 1
    m_Pl = sqrt(8*pi); % Planck mass
    M_Pl = 1;  % Reduced Planck mass
    
end

methods (Static)
    
    function [datastruct] = read_output_file(outfile,record_flags,init_size)
        
        fid = fopen(outfile,'r');
        
        %% Collect meta-data
        
        meta_line = fgets(fid); is = 1;
        
        try
            [n_iter,~,~,is1] = sscanf(meta_line(is:end),'%E,',1); is = is+is1-1;
        catch me
            datastruct = struct();
            return
        end
        
        [mv_0,~,~,is1]   = sscanf(meta_line(is:end),'%G,',1); is = is+is1-1;
        [mh,~,~,is1]     = sscanf(meta_line(is:end),'%G,',1); is = is+is1-1;
        [m_Pl,~,~,is1]    = sscanf(meta_line(is:end),'%G,',1); is = is+is1-1;
        [kmax,~,~,is1]   = sscanf(meta_line(is:end),'%d,',1); is = is+is1-1;
        [gamma,~,~,is1]  = sscanf(meta_line(is:end),'%f,',1); is = is+is1-1;
        
        [measure,~,~,is1]           = sscanf(meta_line(is:end),'%c,',1); is = is+is1-1;
        [n_tunnel_max,~,~,is1]      = sscanf(meta_line(is:end),'%d,',1); is = is+is1-1;
        [lambdascreen,~,~,is1]      = sscanf(meta_line(is:end),'%d,',1); is = is+is1-1;
        [mv_offset_max,~,~,is1]   = sscanf(meta_line(is:end),'%G,',1); is = is+is1-1;
        [fixQ,~,~,is1]              = sscanf(meta_line(is:end),'%d,',1); is = is+is1-1;
        [Nafter,~,~,is1]            = sscanf(meta_line(is:end),'%G,',1); is = is+is1-1;
        [seed,~,~,is1]              = sscanf(meta_line(is:end),'%d,',1); is = is+is1-1;
        [n_recycle,~,~,is1]         = sscanf(meta_line(is:end),'%d,',1); is = is+is1-1;
        
        %% Read output from simulations
        
        if nargin < 3, init_size = round(abs(n_iter)/10); end
        
        data = nan(init_size,16);
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
                    data(i_ln,1:4)  = fscanf(fid,'%G,%d,%G,%G\r\n',4);
                case 2
                    data(i_ln,1:6)  = fscanf(fid,'%G,%d,%G,%G,%d,%G\r\n',5);
                case 3
                    data(i_ln,1:16) = fscanf(fid,'%G,%d,%G,%G,%d,%G,%d,%G,%d,%G,%G,%G,%G,%G,%G,%G\r\n',16);
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
            'mv_offset_max',  mv_offset_max,...
            'fixQ',             logical(fixQ),...
            'Nafter',           Nafter,...
            'seed',             seed,...
            'n_recycle',        n_recycle,...
            'data',             data);
        
        fclose(fid);
        
    end
    
    function [datastruct] = read_output_file_old(outfile,record_flags,init_size)
        
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
        [mv_offset_max,~,~,is1]   = sscanf(meta_line(is:end),'%G,',1); is = is+is1-1;
        [fixQ,~,~,is1]              = sscanf(meta_line(is:end),'%d,',1); is = is+is1-1;
        [Nafter,~,~,is1]            = sscanf(meta_line(is:end),'%G,',1); is = is+is1-1;
        [seed,~,~,is1]              = sscanf(meta_line(is:end),'%d,',1); is = is+is1-1;
        [n_recycle,~,~,is1]         = sscanf(meta_line(is:end),'%d,',1); is = is+is1-1;
        
        %% Read output from simulations
        
        if nargin < 3, init_size = round(abs(n_iter)/10); end
        
        data = nan(init_size,16);
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
                    data(i_ln,1:15) = fscanf(fid,'%G,%d,%G,%d,%G,%d,%G,%d,%G,%G,%G,%G,%G,%G,%G\r\n',15);
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
            'mv_offset_max',  mv_offset_max,...
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

end

% % % Define potential with rescaled mv = 1
% % [V_rescale,Vp_rescale,Vpp_rescale] = build_potential(f,1,mh*obj.m_Pl);
% % f_offset = rho_offset/mv^4/obj.m_Pl^4;
% % V_rescale = @(x) V_rescale(x) - f_offset;
% %
% % B(lr) = fvi.find_tunneling_suppression(R,Y)/mv^4/obj.m_Pl^4;

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
        
