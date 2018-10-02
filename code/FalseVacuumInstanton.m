classdef FalseVacuumInstanton

%% Properties, constructor, and set/get methods

properties ( SetAccess = immutable )
    
    % Potential function and derivatives
    V   = function_handle.empty(); % V(\phi)
    dV  = function_handle.empty(); % \frac{dV}{d\phi}
    d2V = function_handle.empty(); % \frac{d^2V}{d\phi^2}
    
    phi_metaMin         % Location of meta-stable minimum (initial)
    phi_absMin          % Location of absolute minimum    (final)
    
    alpha   = 3         % Dimension of volume-bounding n-sphere (d-1)
    kappa   = 8*pi      % $\kappa = 8\pi G = 8\pi/m_Pl^2 = 1/M_Pl$; Default m_Pl = 1
    phi_bar             % Edge of barrier (where V(phi_bar) = V(phi_metaMin))
    phi_bar_top         % Top of the barrier
    rscale              % Characteristic length scale
    
    no_gravity = false  % Run without the effects of gravity
    flag_plot  = false  % Plot results for debugging
    
end

methods
    
    function self = FalseVacuumInstanton(varargin)
        
        if isscalar(varargin)
            val = varargin{1};
            if ~isstruct(val)
                error('Scalar argument to setParameters must be a struct');
            end
        else
            val = struct(varargin{:});
        end
        
        %% Validate parameters
        
        if isa(val.V,'function_handle')
            self.V = val.V;
        else
            error('V must be a function handle.');
        end
        
        if isa(val.dV,'function_handle')
            self.dV = val.dV;
        else
            error('dV must be a function handle.');
        end
        
        if isa(val.d2V,'function_handle')
            self.d2V = val.d2V;
        else
            error('d2V must be a function handle.');
        end
        
        if isscalar(val.phi_absMin) && isreal(val.phi_absMin)
            self.phi_absMin = val.phi_absMin;
        else
            error('phi_absMin must be a real number.');
        end
        
        if isscalar(val.phi_metaMin) && isreal(val.phi_metaMin)
            self.phi_metaMin = val.phi_metaMin;
        else
            error('phi_metaMin must be a real number.');
        end
        
        % Check that meta-stable vacuum is meta-stable
        if self.V(self.phi_metaMin) <= self.V(self.phi_absMin)
            error('FalseVacuumInstanton:FalseVacuumIsTrue',['V(phi_metaMin) <= V(phi_absMin); tunneling cannot occur.']);
        end
        
        if isfield(val,'alpha')
            if isscalar(val.alpha) && val.alpha > 0 && mod(val.alpha,1) == 0
                self.alpha = val.alpha;
            else
                error('alpha must be a positive integer.');
            end
        end
        
        if isfield(val,'m_Pl')
            if isscalar(val.m_Pl) && isreal(val.m_Pl) && val.m_Pl > 0
                self.kappa = 8*pi/val.m_Pl^2;
            else
                error('m_Pl must be a positive real number');
            end
        end
        
        if isfield(val,'no_gravity')
            if islogical(val.no_gravity) && isscalar(val.no_gravity)
                self.no_gravity = val.no_gravity;
            else
                error('no_gravity must be a scalar logical.');
            end
        end
        
        if isfield(val,'flag_plot')
            if islogical(val.flag_plot) && isscalar(val.flag_plot)
                self.flag_plot = val.flag_plot;
            else
                error('flag_plot must be a scalar logical.');
            end
        end
        
        %% Find barrier location
        
        phi_tol = abs(self.phi_metaMin - self.phi_absMin) * 1e-6;
        V_phimeta = self.V(self.phi_metaMin);
        
        phi1 = self.phi_metaMin;
        phi2 = self.phi_absMin;
        
        phi0 = 0.5 * (phi1+phi2); % Initial guess
        % Do a very simple binary search to narrow down on the right answer.
        while abs(phi1-phi2) > phi_tol
            V0 = self.V(phi0);
            if (V0-V_phimeta) > 0
                phi1 = phi0;
            else
                phi2 = phi0;
            end
            phi0 = 0.5 * (phi1+phi2);
        end
        
        self.phi_bar = phi0;
        
        %% Find the top of the barrier
        
        if isfield(val,'phi_bar_top') 
            if isscalar(val.phi_bar_top) && isreal(val.phi_bar_top)
                self.phi_bar_top = val.phi_bar_top;
                phi_top = self.phi_bar_top;
            else
                error('phi_bar_top must be a real number.');
            end
        else
            phi_tol = abs(self.phi_bar - self.phi_metaMin) * 1e-6;
            phi1 = self.phi_metaMin;
            phi2 = self.phi_bar;
            sgn = sign(self.dV(self.phi_bar));
            phi_top = 0.5 * (phi1+phi2);
            while abs(phi1-phi2) > phi_tol
                if sgn*self.dV(phi_top) < 0
                    phi1 = phi_top;
                else
                    phi2 = phi_top;
                end
                phi_top = 0.5 * (phi1 + phi2);
            end
            self.phi_bar_top = phi_top;
        end
        
        if ~((self.phi_bar < phi_top && phi_top < self.phi_metaMin) || ...
                (self.phi_bar > phi_top && phi_top > self.phi_metaMin))
            error(['Minimization is placing the top of the ' ...
                'potential barrier outside of the interval defined by ' ...
                'phi_bar and phi_metaMin. Assume that the barrier does not exist.' ...
                'no barrier']);
        end
        
        %% Find characteristic length scale
        
        Vtop = self.V(phi_top) - self.V(self.phi_metaMin);
        xtop = phi_top - self.phi_metaMin;
        % Cubic function given by (ignoring linear and constant terms):
        % f(x) = C [(-1/3)x^3 + (1/2)x^2 xtop]
        % C = 6 Vtop / xtop^3
        % f''(xtop) = - C xtop
        % d2V = -6*Vtop / xtop^2
        % rscale = 1 / sqrt(d2V)
        if Vtop <= 0
            error(['Barrier height is not positive, ' ...
                'does not exist.', 'no barrier']);
        end
        self.rscale = abs(xtop) / sqrt(abs(6*Vtop));
        
    end
    
end

%% Integration

methods (Access = protected)
    
    function [r0,phi_r0,dphi_r0] = initial_conditions_no_gravity(self,delta_phi0,rmin,delta_phi_cutoff)
        % Finds the initial conditions for integration.
        % 
        % The instanton equations of motion are singular at `r=0`, so we
        % need to start the integration at some larger radius. This
        % function finds the value `r0` such that `phi(r0) = phi_cutoff`.
        % If there is no such value, it returns the initial conditions at `rmin`.
        
        %% Try out guess given by delta_phi0
        
        phi0 = self.phi_absMin + delta_phi0;
        dV0  = self.dV(phi0);
        d2V0 = self.d2V(phi0);
        
        [phi_r0, dphi_r0] = self.exact_solution(rmin, phi0, dV0, d2V0);
        if abs(phi_r0 - self.phi_absMin) > abs(delta_phi_cutoff)
            % The initial conditions at rmin work. Stop here.
            r0 = rmin;
            return
        end
        if sign(dphi_r0) ~= sign(delta_phi0)
            % The field is evolving in the wrong direction.
            % Increasing r0 won't increase |delta_phi_r0|/
            r0 = rmin;
            return
        end
        
        %% Find a value of r0 such that delta_phi_r0 > delta_phi_cutoff
        
        r = rmin;
        while isfinite(r)
            rlast = r;
            r = r*10;
            phi = self.exact_solution(r, phi0, dV0, d2V0);
            if abs(phi - self.phi_absMin) > abs(delta_phi_cutoff)
                break
            end
        end
        
        %% Now find where (phi - self.phi_absMin) = delta_phi_cutoff exactly
        
        deltaPhiDiff = @(r_) feval(@(p) abs(p - self.phi_absMin) - abs(delta_phi_cutoff),...
            self.exact_solution(r_, phi0, dV0, d2V0));
        r0 = fzero(deltaPhiDiff, [rlast, r]);
        [phi_r0, dphi_r0] = self.exact_solution(r0, phi0, dV0, d2V0);
        
    end
    
    function [r0,y0] = initial_conditions(self,delta_phi0,rmin,delta_phi_cutoff)
        % Finds the initial conditions for integration.
        % 
        % When including gravity, we also need to find the initial
        % conditions for the angular radius `rho(r0)`. However, the coupled
        % set of differential equations for `rho` and `phi` is non-linear,
        % so we cannot plausibly find an analytic solution. Instead, this
        % function assumes that `dphi/dr` is small and `V(phi)` is
        % approximately constant for the purpose of calculating `rho(r)`,
        % which can then be found analytically. An initial guess `r0_guess`
        % is calculated using
        % :method:`SingleFieldInstanton.initial_conditions`. A new initial
        % radius is then chosen such that `(1/rho(r0))(drho/dr) =
        % 1/r0_guess`. This way, the value of the friction term at `r0` is
        % the same as it would be without gravity.
        
        %% Get initial guess without gravity
        
        [r0_guess,phi0,dphi0] = initial_conditions_no_gravity(...
            self,delta_phi0,rmin,delta_phi_cutoff);
        
        V0 = self.V(phi0);
        
        %% Choose rho0, drho0 s.t. r0_guess = rho/drho
        
        w = abs(self.kappa/3*V0)^0.5;
        if w == 0 || self.no_gravity
            % No curvature.
            r0    = r0_guess;
            rho0  = r0_guess;
            drho0 = 1;
        elseif V0 > 0
            % positive curvature
            r0    = atan(w*r0_guess)/w;
            rho0  = sin(w*r0)/w;
            drho0 = cos(w*r0);
        else
            % negative curvature
            if w*r0_guess >= 1
                error('FalseVacuumInstanton:StableFalseVacuum',[...
                    'The friction term in the negatively curved space never ' ...
                    'gets small enough to allow tunneling to proceed.'...
                    'stable, not metastable']);
            end
            r0    = atanh(w*r0_guess)/w;
            rho0  = sinh(w*r0)/w;
            drho0 = cosh(w*r0);
        end
        
        y0 = [phi0,dphi0,rho0,drho0];
        
    end
    
    function [r1,y1,convergence_type] = integrate_profile(self,r0,y0,epsabs,rmax)
        % Integrate the bubble wall equation
        % 
        % This works the same basic way as when there is no gravity, except
        % now the overshoot/undershoot conditions are a little bit more
        % complicated. The solution will overshoot if at any point `rho(r)
        % = 0` for `r > 0` (indicating that the bubble has wrapped around
        % to the anti-podal point of de Sitter space), unless `dphi/dr` is
        % also zero.
        
        dr = rmax;
        
        % dY is the ODE that we use
        dY = @(r,y) self.equation_of_motion(r,y).';
        
        % Positive means we're heading down
        ysign = sign(y0(1)-self.phi_metaMin);
        
        % Perform integration using stiff ODE solver
        options = odeset(...
            'Jacobian', @(t,y) self.eom_jacobian(t,y),...
            'Events',   @(~,y) self.ode_events(y,epsabs,ysign) );
        [r1,y1,re,~,ie] = ode23s(dY,[r0,r0+dr],y0.',options);
        
        if ~any(re == r1(end))
            error('FalseVacuumInstanton:NoTerminatingEvents',...
                'No terminating events were triggered.')
        end
        
        switch ie(end)
            case {1,2,6} % Got close enough to meta-stable minimum with small enough derivative
                convergence_type = 'converged';
            case {3} % Didn't make it over the potential barrier
                convergence_type = 'undershoot';
            case {4} % Passed meta-stable minimum without stopping
                convergence_type = 'overshoot';
            case {5}
                if abs(y1(end,2)) < epsabs(2)
                    convergence_type = 'converged';
                else
                    convergence_type = 'overshoot';
                end
            otherwise
                error('FalseVacuumInstanton:IntegralDiverged',...
                    'ODE solver failed to integrate instanton solution.');
        end
        
        if self.flag_plot
            w = abs(self.kappa/3*self.V(self.phi_absMin))^0.5;
            subplot(2,2,1); plot(r1*w,y1(:,1));
            subplot(2,2,2); plot(r1*w,y1(:,2));
            subplot(2,2,3); plot(r1*w,y1(:,3));
            subplot(2,2,4);
            x = linspace(y1(end,1),y1(1,1),1000); plot(x,self.V(x)); hold on;
            plot(y1(end,1),self.V(y1(end,1)),'o');
            plot(y1(1,1),self.V(y1(1,1)),'xk'); hold off
            disp(num2str(ie.')); pause(0.5); 
        end
        
    end
    
end

%% Physics

methods
    
    function dy = equation_of_motion(self,~,y)
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
        rho  = y(:,3);
        drho = y(:,4);
        
        % EOM for the inflaton field
        d2phi = self.dV(phi.').' - self.alpha*drho.*dphi./rho;
        
        % EOM for spacetime geometry
        if self.no_gravity
            d2rho = zeros(size(drho));
        else
            d2rho = -self.kappa/3*rho.*(dphi.^2 + self.V(phi.').');
        end
        
        dy = [dphi,d2phi,drho,d2rho];
        
    end
    
    function jac = eom_jacobian(self,~,y)
        
        % Reshape the input if necessary
        y = permute(y,circshift([1 2],[0,find(size(y) == 4)]));
        
        % Unpack variables
        phi  = y(:,1);
        dphi = y(:,2);
        rho  = y(:,3);
        drho = y(:,4);
        
        jac = zeros(4);
        
        % Assemble Jacobian
        if ~self.no_gravity
            jac(1,:) = [0 1 0 0];
            jac(2,:) = [self.d2V(phi), -self.alpha*drho./rho, self.alpha*drho.*dphi./rho^2, -self.alpha*dphi./rho];
            jac(3,:) = [0 0 0 1];
            jac(4,:) = [rho*self.dV(phi), 2*rho*dphi, (dphi.*dphi + self.V(phi)), 0]*(-self.kappa/3);
        else
            jac(1,:) = [0 1 0 0];
            jac(2,:) = [self.d2V(phi), -self.alpha*drho./rho, 0, 0];
            jac(3,:) = [0 0 0 1];
            jac(4,:) = [0 0 0 0];
        end
        
    end
    
    function [phi,dphi] = exact_solution(self,r,phi0,dV,d2V)
        % Find `phi(r)` given `phi(r=0)`, assuming a quadratic potential.
        % 
        % If the potential at the point :math:`\phi_0` is a simple quadratic, the
        % solution to the instanton equation of motion can be determined exactly.
        % The non-singular solution to 
        % 
        % .. math::
        %   \frac{d^2\phi}{dr^2} + \frac{\alpha}{r}\frac{d\phi}{dr} =
        %   V'(\phi_0) + V''(\phi_0) (\phi-\phi_0)
        % 
        % is
        % 
        % .. math::
        %   \phi(r)-\phi_0 = \frac{V'}{V''}\left[
        %   \Gamma(\nu+1)\left(\frac{\beta r}{2}\right)^{-\nu} I_\nu(\beta r) - 1
        %   \right]
        % 
        % where :math:`\nu = \frac{\alpha-1}{2}`, :math:`I_\nu` is the modified
        % Bessel function, and :math:`\beta^2 = V''(\phi_0) > 0`. If instead 
        % :math:`-\beta^2 = V''(\phi_0) < 0`, the solution is the same but with 
        % :math:`I_\nu \rightarrow J_\nu`.
        
        beta   = sqrt(abs(d2V));
        beta_r = beta*r;
        nu     = 0.5 * (self.alpha - 1);
        
        if d2V > 0
            bessfun = @besseli;
            s = +1;
        else
            bessfun = @besselj;
            s = -1; 
        end
        
        phi  = dV/d2V*(gamma(nu+1)*(0.5*beta_r)^(-nu)*bessfun(nu,beta_r)-1) + phi0;
        
        dphi = dV/d2V*gamma(nu+1)*(0.5*beta_r)^(-nu)*(...
            (-nu/r) * bessfun(nu,beta_r) + ...
            0.5*beta * (bessfun(nu-1,beta_r)+s*bessfun(nu+1,beta_r)) );
        
    end
    
    function [value,isterminal,direction] = ode_events(self,y,eps,ysign)
        
        phi = y(1); dphi = y(2);
        phi0 = self.phi_metaMin;
        
        if self.no_gravity
            
            % 1  Possibly converged; phi = phi_metaMin within tolerance
            % 2  Possibly converged; dphi = 0 within tolerance
            %       Terminate only if both 1 & 2 are satisfied
            % 3  Undershoot; dphi/dr has changed sign, not converged
            % 4  Overshoot; gone past phi_metaMin beyond tolerance
            % 5  Overshoot; dphi/dr has likely diverged
            
            events = [...
                (phi-phi0)*ysign - eps(1),  abs(dphi) < eps(2),         -1;...
                dphi*ysign + eps(2),        abs(phi-phi0) < eps(1),     +1;...
                dphi*ysign - eps(2),        (phi-phi0)*ysign > eps(1),  +1;...
                (phi-phi0)*ysign + eps(1),  true,                       -1
            ];
            
        else
            
            % 1  Possibly converged; drho/dr close to -1, so rho close to 0
            %       (Easier to converge than demanding rho = 0)
            %       Terminate iff dphi = 0
            % 2  Possibly converged; phi = phi_metaMin within tolerance
            %       Terminate iff dphi = 0
            % 3  Possibly converged; dphi = 0 within tolerance
            %       Terminate iff phi = phi_metaMin
            % 4  Undershoot; dphi/dr has changed sign, not converged
            % 5  Overshoot; gone past phi_metaMin beyond tolerance
            % 6  Overshoot; Collapsed on rho = 0, but abs(dphi/dr) > 0
            %       Equations of motion are singular, not a valid solution
            
            rho = y(3); drho = y(4);
            
            events = [...
                (phi-phi0)*ysign - eps(1),  abs(dphi) < eps(2),         -1;...
                dphi*ysign + eps(2),        (phi-phi0)*ysign < eps(1) || ...
                                                (drho < -1+eps(4)),     +1;...
                dphi*ysign - eps(2),        (phi-phi0)*ysign > eps(1)   +1;...
                (phi-phi0)*ysign + eps(1),  true,                       -1;...
                drho - (-1+eps(4)),         true                        -1;... % abs(dphi) < eps(2)
                rho,                        true,                       -1;...
                ...%abs(dphi/phi) - 1e5,        true,                       +1
            ];
            
        end
        
        value       = events(:,1);
        isterminal  = events(:,2);
        direction   = events(:,3);
        
    end
    
end

%% Main

methods
    
    function [R,Y,useThinWall] = find_profile(self,xguess,xtol,phitol,thinCutoff,B_cutoff,rmin,rmax)
        % Calculate the bubble profile by iteratively over/undershooting.
        % 
        % This is very similar to :method:`SingleFieldInstanton.find_profile`,
        % but slightly modified for the inclusion of gravity.
        %
        % This will call :func:`integrate_profile` many times, trying to find
        % the correct initial condition `phi(r=0)` such that the field ends up
        % in the metastable vacuum at infinity. Once the correct initial
        % condition is found, it calls :func:`integrate_and_save_profile` to find
        % the profile along the length of the wall.
        % 
        % For very thin-walled bubbles, the initially value of `phi` can be
        % extremely close to the stable minimum and small variations in `phi`
        % can cause large variations in the integration. Rather than varying 
        % `phi(r=0)` directly, it is easier to vary a parameter `x` defined by
        % 
        % .. math::
        %    \phi(r=0) = \phi_{\rm absMin} 
        %    + e^{-x}(\phi_{\rm metaMin}-\phi_{\rm absMin})
        % 
        % This way, `phi = phi_metaMin` when `x` is zero and 
        % `phi = phi_absMin` when `x` is  infinity.
        
        %% Set bounds and error tolerances
        
        if nargin < 2, xguess           = [];   end
        if nargin < 3, xtol             = 1e-4; end
        if nargin < 4, phitol           = 1e-4; end
        if nargin < 5, thinCutoff       = 1e-2; end
        if nargin < 6, B_cutoff         = Inf;  end
        if nargin < 6, rmin             = 1e-4; end
        if nargin < 7, rmax             = 1e+4; end
        
        delta_phi = self.phi_metaMin - self.phi_absMin;
        
        % Set x parameters
        xtop = -log(abs((self.phi_bar_top-self.phi_absMin)/delta_phi));
        xbar = -log(abs((self.phi_bar-self.phi_absMin)/delta_phi));
        
        xmin = xtop; % No further toward meta min than top of the barrier
        xmax = Inf;  % No initial limit on how close to absolute minimum
        
        if ~isempty(xguess)
            x = xguess;
        else
            x = xbar;
        end
        xincrease = 2.0;
        
        % Set r parameters
        rmin  = rmin * self.rscale;
        rmax  = rmax * self.rscale;
        
        % Set the phi parameters
        delta_phi_cutoff = thinCutoff * delta_phi;
        epsabs = abs([...
            delta_phi*phitol,...
            delta_phi*phitol/self.rscale,...
            0.1 * rmin,...
            0.005 ]);
        
        %% Repeat profile integration until it converges
        
        ctype_last = ''; dphi_last = Inf;
        
        useThinWall = false;
        
        maxloops = 100;
        dysign = nan;
        R = NaN;
        for nloops = 1:maxloops
            
            % Set initial conditions
            delta_phi0 = exp(-x)*delta_phi;
            [r0,y0] = self.initial_conditions(delta_phi0,rmin,delta_phi_cutoff);
            if ~isfinite(r0) || ~isfinite(x)
                % Use the last finite values instead (assuming there are such values)
                assert(~isnan(R),'Failed to retrieve initial conditions on the first try.');
                break
            end
            
            if isnan(dysign)
                dysign = sign(y0(2));
            elseif sign(y0(2)) ~= dysign
                % Initialized going the wrong way. The true location of the
                % minimum is slightly off from phi_absMin, and we've
                % crossed it. Treat as overshoot and try again.
                xmax = x;
                x = .5*(xmin+xmax);
                continue
            end
            
            % Integrate profile
            [R,Y,ctype] = self.integrate_profile(r0,y0,epsabs,rmax);
            
            %% Check for Hawking-Moss instanton
            
%             if xmin == xbar && ...
%                     strcmp(ctype_last,'overshoot') && ...
%                     abs(dphi_last) < abs(Y(end,2))-epsabs(2)
%                 % Passed minimum of abs(dphi) \rvert_{phi = phi_metaMin}
%                 % There is no CDL bubble; the only solution is the Hawking-
%                 % Moss instanton that sits on top of the barrier forever.
%                 w_top = abs(self.kappa/3*self.V(self.phi_bar_top))^0.5;
%                 R = pi/w_top/2;
%                 Y = [self.phi_bar_top,0,1/w_top,0,-w_top];
%                 useThinWall = false;
%                 return
%             end
            
            ctype_last = ctype; dphi_last = Y(end,2);
            
            %% Check for overshoot, undershoot
            
            switch ctype
                case'converged'
                    break
                case 'undershoot' % x is too low
                    xmin = x;
                    if isinf(xmax)
                        x = x*xincrease;
                    else
                        x = .5*(xmin+xmax);
                    end
%                     xx = linspace(0e3,2e3,1001);
%                     plot(xx,self.V(xx)); hold on;
%                     plot(exp(-x)*delta_phi+self.phi_absMin,self.V(exp(-x)*delta_phi+self.phi_absMin),'o'); hold off;
%                     pause(0.5);
                case 'overshoot' % x is too high
                    xmax = x;
                    x = .5*(xmin+xmax);
            end
            
            % Close enough; don't wait for convergence
            if (xmax-xmin) < xtol
                if abs(xmax-xtop) < xtol
                    % Passed minimum of abs(dphi) \rvert_{phi = phi_metaMin}
                    % There is no CDL bubble; the only solution is the Hawking-
                    % Moss instanton that sits on top of the barrier forever.
                    w_top = abs(self.kappa/3*self.V(self.phi_bar_top))^0.5;
                    R = pi/w_top/2;
                    Y = [self.phi_bar_top,0,1/w_top,0,-w_top];
                    useThinWall = false;
                    return
                end
                break
            end
            
%             if ~isinf(self.B_cutoff)
%                 B = self.find_tunneling_suppression(R,Y);
%                 disp(['B = ' num2str(B)]);
%                 if B > B_cutoff
%                     Y(1,1) = nan;
%                     return
%                 end
%             end
             
        end
        
        % Add d^2Rho/dr^2 to the profile
        dYdR = self.equation_of_motion(R,Y);
        Y = [Y dYdR(:,end)];
        
        return
        
        % Recommend thin wall approximation if the bubble size is much less
        % than the de Sitter radius of the background.
        w = abs(self.kappa/3*self.V(self.phi_metaMin))^0.5;
        if w*R(end) < 1e-2, useThinWall = true; end
        
    end
    
    function S = find_action_4D(self,R,Y,actionForm)
        
        if self.no_gravity
            actionForm = 0;
        elseif nargin < 4
            actionForm = 2;
        end
        
        S = self.find_action_4D_static(self.alpha,self.kappa,self.V,R,Y,actionForm);
        
    end
    
    function B = find_tunneling_suppression(self,R,Y,actionForm)
        % Computes the exponential suppression of the tunneling rate
        %
        % If \Gamma/V = A*exp(-B/\hbar)*(1 + O(\hbar)), then B is the
        % exponential suppression of the tunneling rate
        
        if nargin < 4
            B = find_tunneling_suppression_static(self.alpha,self.kappa,...
                self.V,self.phi_metaMin,R,Y);
        else
            if self.no_gravity
                actionForm = 0;
            end
            B = find_tunneling_suppression_static(self.alpha,self.kappa,...
                self.V,self.phi_metaMin,R,Y,actionForm);
        end
        
    end
    
    function B = parke_thinwall(self,R,Y)
        
        Vf = self.V(self.phi_metaMin);
        Vt = self.V(self.phi_absMin);
        
        if nargin < 3
            % Assume the form of the profile when calculating bubble tension
            phi = @(r) (self.phi_metaMin - self.phi_absMin)/2*(1+tanh(r/self.rscale/1.8)) + self.phi_absMin;
            integrand = @(r) self.V(phi(r)) - Vf + (phi(r) - self.phi_metaMin)*(Vf-Vt)/(self.phi_absMin - self.phi_metaMin);
            S1 = 2*integral(integrand,-4*self.rscale,4*self.rscale);
        else
            % Use the supplied profile
            Y = mat2cell(Y,size(Y,1),ones(1,size(Y,2))); Phi = Y{1};
            S1 = 2*trapz(R, self.V(Phi) - Vf + (Phi - self.phi_metaMin)*...
                (Vf - Vt)/(self.phi_absMin - self.phi_metaMin));
        end
        
        rhobar0 = 3*S1/(Vf-Vt);
        lambda  = (self.kappa/3*(Vf+Vt))^(-1/2);
        Lambda  = (self.kappa/3*(Vf-Vt))^(-1/2);
        
        % Angular radius at bubble wall (approx. = R0)
        rhobar  = rhobar0/sqrt(1 + 2*(rhobar0/2/lambda)^2 + (rhobar0/2/Lambda)^4);
        
        % Check self-consistency of thin wall approximation
        thin_thres = 0.01;
        if nargin < 3 && self.rscale/rhobar > thin_thres
            B = nan;
            return
        end
        
        % Contribution from wall
        B_wall = 2*pi^2*rhobar^3*S1;
        
        % Contribution from bubble interior
        B_inside = 12*pi^2/self.kappa*(...
            ((1-self.kappa/3*Vt*rhobar^2)^(3/2)-1)/Vt - ...
            ((1-self.kappa/3*Vf*rhobar^2)^(3/2)-1)/Vf);
        
        B = B_inside + B_wall;
        
        return
        
        %% Compute pre-factor
        
        % Calculate pre-factor sans Hubble factor exp(3Ht)
        % exp(zeta_R'(-2)) = exp(-zeta_R(3)/(2*pi)^2) = 0.97001
        B = B - log(0.25*0.97001*rhobar^2*S1^2);
        % are the units of rhobar and S1 correct?
        
    end
    
    function B = parke_thinwall2(self)
        
        Vf = self.V(self.phi_metaMin);
        Vt = self.V(self.phi_absMin);
        
        % Assume the form of phi when calculating bubble tension
        phi = @(r) (self.phi_metaMin - self.phi_absMin)/2*(1+tanh(r/self.rscale)) + self.phi_absMin;
        integrand = @(r) self.V(phi(r)) - Vf + (phi(r) - self.phi_metaMin)*(Vf-Vt)/(self.phi_absMin - self.phi_metaMin);
        S1 = 2*integral(integrand,-4*self.rscale,4*self.rscale);
        
        rhobar0 = 3*S1/(Vf-Vt);
        lambda  = (self.kappa/3*(Vf+Vt))^(-1/2);
        Lambda  = (self.kappa/3*(Vf-Vt))^(-1/2);
        
        % Angular radius at bubble wall (approx. = R0)
        rhobar  = rhobar0/sqrt(1 + 2*(rhobar0/2/lambda)^2 + (rhobar0/2/Lambda)^4);
        
        % Contribution from wall
        B_wall = 2*pi^2*rhobar^3*S1;
        
        % Contribution from bubble interior
        B_inside = 12*pi^2/self.kappa*(...
            ((1-self.kappa/3*Vt*rhobar^2)^(3/2)-1)/Vt - ...
            ((1-self.kappa/3*Vf*rhobar^2)^(3/2)-1)/Vf);
        
        B = B_inside + B_wall;
        
        %% Compute pre-factor
        
        % Calculate pre-factor sans Hubble factor exp(3Ht)
        % exp(zeta_R'(-2)) = exp(-zeta_R(3)/(2*pi)^2) = 0.97001
%         B = B - log(0.25*0.97001*rhobar^2*S1^2);
        % are the units of rhobar and S1 correct?
        
    end
    
end

methods (Static)
    
    function B = find_tunneling_suppression_static(alpha,kappa,V,phi_metaMin,R,Y,actionForm)
        % Computes the exponential suppression of the tunneling rate
        %
        % If \Gamma/V = A*exp(-B/\hbar)*(1 + O(\hbar)), then B is the
        % exponential suppression of the tunneling rate
        
        if nargin < 7, actionForm = 2; end
        
        Y = mat2cell(Y,size(Y,1),ones(1,size(Y,2)));
        
        %% Set up queries for the false vacuum background
        % The "interior" and "exterior" volume terms are all we need, so
        % just provide a single data point specifying phi_metaMin.
        % find_action_4D won't bother integrating the Lagrangian.
        
        w = abs(kappa/3*V(phi_metaMin))^0.5;
        if w == 0 || no_gravity
            % Integrate inside the full support of the bubble.
            R_bkgd     = R(end);
            Rho_bkgd   = R(end);
            dRho_bkgd  = 1;
            d2Rho_bkgd = 0;
        elseif V(phi_metaMin) > 0
            % Integrate over the whole compact space. (Any R would work.)
            R_bkgd     = pi/2/w;
            Rho_bkgd   = 1/w;
            dRho_bkgd  = 0;
            d2Rho_bkgd = -w;
        else
            % Integrate until the radius of curvature of the background
            % space matches that at the outer edge of the bubble wall.
            R_bkgd     = asinh(Y{3}(end)*w)/w;
            Rho_bkgd   = Y{3}(end);
            dRho_bkgd  = cosh(R_bkgd*w);
            d2Rho_bkgd = sinh(R_bkgd*w)*w;
        end
        
        Y_bkgd = {phi_metaMin,0,Rho_bkgd,dRho_bkgd,d2Rho_bkgd};
        
        %% Compute tunneling rate exponential piece
        
        S_bubble = FalseVacuumInstanton.find_action_4D_static(alpha,kappa,V,R,Y,actionForm);
        S_bkgd   = FalseVacuumInstanton.find_action_4D_static(alpha,kappa,V,R_bkgd,Y_bkgd,actionForm);
        
        B = S_bubble-S_bkgd;
        
    end
    
    function S = find_action_4D_static(alpha,kappa,V,R,Y,actionForm)
        
        if nargin < 4
            actionForm = 2;
        end
        
        [Phi,dPhi,Rho,dRho,d2Rho] = deal(Y{:});
        
        %% Compute volume terms for bubble interior/exterior (assumes 4D)
        % Generalize to d dimensions?
        
        w_int = abs(kappa/3*V(Phi(1)))^0.5;
        
        % Integrate over the bubble interior
        if w_int == 0 || actionForm == 0
            
            % field part    : integrate 2*pi^2*r^3*V dr from 0 to R
            volume = pi^2*R(1)^4/2;
            S = volume*V(Phi(1));
            
        elseif V(Phi(1)) > 0
            
            % field part    : integrate 2*pi^2*sin(w*r)^3/w^3*V dr from 0 to R
            % geometry part : integrate 2*pi^2/k*(-sin(w*r)^3/w + sin(w*r)*cos(w*r)^2/w - sin(w*r)/w) dr from 0 to R
            volume = 8*pi^2/3/w_int^4*sin(R(1)*w_int/2)^4*(cos(R(1)*w_int)+2); % = pi^2/6/w^4*(cos(3*R(1)*w)-9*cos(R(1)*w));
            S = volume*(V(Phi(1)) - 6*w_int^2/kappa);
            
        else
            
            % field part    : integrate 2*pi^2*sinh(w*r)^3/w^3*V dr from 0 to R
            % geometry part : integrate 2*pi^2/k*(sinh(w*r)^3/w + sinh(w*r)*cosh(w*r)^2/w - sinh(w*r)/w) dr from 0 to R
            volume = 8*pi^2/3/w_int^4*sinh(R(1)*w_int/2)^4*(cosh(R(1)*w_int)+2); % = pi^2/6/w^4*(cosh(3*R(1)*w)-9*cosh(R(1)*w));
            S = volume*(V(Phi(1)) + 6*w_int^2/kappa);
            
        end
        
        if V(Phi(end)) > 0 || actionForm ~= 0
            
            % Integrate over the finite volume outside the bubble
            w_ext = abs(kappa/3*V(Phi(end)))^0.5;
            Rend = asin(w_ext*Rho(end))/w_ext;
            if dRho(end) < 0, Rend = pi/w_ext - Rend; end
            volume = 8*pi^2/w_ext^4 - 8*pi^2/3/w_ext^4*sin(Rend*w_ext/2)^4*(cos(Rend*w_ext)+2);
            S = S + volume*(V(Phi(end)) - 6*w_ext^2/kappa);
            
        end
        
        %% Integrate Lagrangian for bubble wall
        
        if isscalar(R), return, end
        
        % Find the area of an n-sphere (alpha=n):
        d = alpha + 1; % Number of dimensions in the integration
        area = Rho.^alpha * 2*pi^(d*.5)/gamma(d*.5);
        
        % Select Lagrangian
        switch actionForm
            case 0
                % Matter-only Lagrangian -- no gravity
                lagr = 0.5*dPhi.^2 + V(Phi);
                boundary = 0;
            case 1
                % Full Einstein-Hilbert Lagrangian
                lagr = 0.5*dPhi.^2 + V(Phi) + 3./kappa.*(d2Rho./Rho + (dRho./Rho).^2 - Rho.^(-2));
                boundary = 0;
            case 2
                % After integration by parts
                lagr = 0.5*dPhi.^2 + V(Phi) - 3./kappa.*((dRho./Rho).^2 + Rho.^(-2));
                boundary = 6*pi^2/kappa*(Rho(end)^2*dRho(end) - Rho(1)^2*dRho(1));
            case 3
                % After integration by parts and asserting on-shell GR
                lagr = 2*V(Phi) - 6./kappa./Rho.^2;
                boundary = 6*pi^2/kappa*(Rho(end)^2*dRho(end) - Rho(1)^2*dRho(1));
        end
        
        %% Compute action for bubble wall
        
        S = S + trapz(R,area.*lagr) + boundary;
        
    end
    
end

end

%% Scratch

% % %% Utilities
% % 
% % methods (Access = protected)
% %     
% %     function dV_phi = dV_near_min(self,delta_phi,phi_min)
% %         % Calculates `dV/dphi` at ``phi = phi_absMin + delta_phi``.
% %         
% %         if nargin < 3, phi_min = self.phi_absMin; end
% %         
% %         phi = phi_min + delta_phi;
% %         dV_phi = self.dV(phi);
% %         
% %     end
% %     
% % end

% % function f = cubicInterpFunction(y0,dy0,y1,dy1)
% %     
% %     y3 = y1;
% %     y1 = y0 + dy0/3.0;
% %     y2 = y3 - dy1/3.0;
% %     
% %     f = @(t) feval(@(mt) y0*mt.^3 + 3*y1*mt.*mt.*t + 3*y2*mt.*t.*t + y3*t.^3,1-t);
% %     
% % end

% %     function S = find_action_without_gravity(self,R,Y)
% %         % Calculate the Euclidean action for the instanton:
% %         % 
% %         % .. math::
% %         %   S = \int [(d\phi/dr)^2 + V(\phi)] r^\alpha dr d\Omega_\alpha
% %         
% %         Y = mat2cell(Y,ones(1,size(Y,1)),size(Y,2));
% %         [Phi,dPhi,~,~,~] = deal(Y{:});
% %         
% %         % Find the area of an n-sphere (alpha=n):
% %         d = self.alpha + 1; % Number of dimensions in the integration
% %         area = R.^self.alpha * 2*pi^(d*.5)/gamma(d*.5);
% %         
% %         % And integrate the profile
% %         lagr = 0.5 * dPhi.^2 + self.V(Phi) - self.V(self.phi_metaMin);
% %         S = trapz(R,lagr.*area);
% %         
% %         % Find the bulk term in the bubble interior
% %         % Assume field is uniform inside the bubble
% %         volume = R(1)^d * pi^(d*.5)/gamma(d*.5 + 1);
% %         S = S + volume * (self.V(Phi(1)) - self.V(self.phi_metaMin));
% %         
% %     end
% %     
% %     function B = find_action_old(self,R,Y)
% %         % Calculate the Euclidean action for the instanton:
% %         % 
% %         % .. math::
% %         %   S = \int [(d\phi/dr)^2 + V(\phi)] r^\alpha dr d\Omega_\alpha
% %         
% %         % f[r_?NumberQ] := - r^(rDim-1) profile[r] . (gradient[profile[r]]-gradient[falseVacuum]);
% %         % Pi^(rDim/2)/Gamma[1+rDim/2] NIntegrate[f[r], Join[{r, 0}, ri, {Infinity}],
% %         
% %         if self.no_gravity
% %             B = self.find_action_without_gravity(R,Y);
% %             return
% %         end
% %         
% %         Y = mat2cell(Y,ones(1,size(Y,1)),size(Y,2));
% %         [Phi,dPhi,Rho,dRho,d2Rho] = deal(Y{:});
% %         
% %         % Find the area of an n-sphere (alpha=n):
% %         d = self.alpha + 1; % Number of dimensions in the integration
% %         area = Rho.^self.alpha * 2*pi^(d*.5)/gamma(d*.5);
% %         
% %         actionForm = 0;
% %         %   0   Full Euclidean EH action
% %         %   1   After integration by parts
% %         %   2   Assuming on-shell GR
% %         
% %         %% Action for bubble
% %         
% %         switch actionForm
% %             case 0
% %                 lagr = 0.5*dPhi.^2 + self.V(Phi) + 1./self.kappa.*(d2Rho./Rho + (dRho./Rho).^2 - Rho.^(-2));
% %             case 1
% %                 lagr = 0.5*dPhi.^2 + self.V(Phi) - 1./self.kappa.*((dRho./Rho).^2 + Rho.^(-2));
% %             case 2
% %                 lagr = 2*self.V(Phi) - 2./self.kappa./Rho.^2;
% %         end
% %         S_bubble = trapz(R,area.*lagr);
% %         
% %         if actionForm > 0
% %             S_bubble = S_bubble + 2*pi^2/self.kappa*(Rho(end)^2*dRho(end));
% %         end
% %         
% %         w = abs(self.kappa*self.V(Phi(1)))^0.5;
% %         
% %         % Add volume term for bubble interior
% %         if self.V(Phi(1)) > 0
% %             % field part    : integrate 2*pi^2*sin(w*r)^3/w^3*V dr from 0 to R
% %             % geometry part : integrate 2*pi^2/k*(-sin(w*r)^3/w + sin(w*r)*cos(w*r)^2/w - sin(w*r)/w) dr from 0 to R
% %             volume = 8*pi^2/3/w^4*sin(R(1)*w/2)^4*(cos(R(1)*w)+2);
% %             S_bubble = S_bubble + volume*(2*self.V(Phi(1))) - ...
% %                 16*pi^2*sin(R(1)*w/2)^4*(cos(R(1)*w)+2)/3/self.kappa/w^2;
% %         else
% %             % field part    : integrate 2*pi^2*sinh(w*r)^3/w^3*V dr from 0 to R
% %             % geometry part : integrate 2*pi^2/k*(sinh(w*r)^3/w + sinh(w*r)*cosh(w*r)^2/w - sinh(w*r)/w) dr from 0 to R
% %         end
% %         
% %         %% Action for false vacuum background
% %         
% %         % Get Rho(R) in the case of uniform false vacuum
% %         w = abs(self.kappa*self.V(self.phi_metaMin))^0.5;
% %         if w == 0
% %             Rho_metaMin = R;
% %         elseif self.V(self.phi_metaMin) > 0
% %             Rho_metaMin = sin(R*w)/w;
% %             dRho_metaMin = cos(R*w);
% %             d2Rho_metaMin = -sin(R*w)*w;
% %         else
% %             Rho_metaMin = sinh(R*w)/w;
% %             dRho_metaMin = cosh(R*w);
% %             d2Rho_metaMin = sinh(R*w)*w;
% %         end
% %         area_metaMin = Rho_metaMin.^self.alpha * 2*pi^(d*.5)/gamma(d*.5);
% %         
% %         switch actionForm
% %             case 0
% %                 lagr = self.V(self.phi_metaMin) + ...
% %                     1./self.kappa.*(d2Rho_metaMin./Rho_metaMin + (dRho_metaMin./Rho_metaMin).^2 - Rho_metaMin.^(-2));
% %             case 1
% %                 lagr = self.V(self.phi_metaMin) - 1./self.kappa.*((dRho_metaMin./Rho_metaMin).^2 + Rho_metaMin.^(-2));
% %             case 2
% %                 lagr = 2*self.V(self.phi_metaMin) - 2./self.kappa./Rho_metaMin.^2;
% %         end
% %         S_metaMin = trapz(R,area_metaMin.*lagr);
% %         
% %         % Add the surface term
% %         if actionForm > 0
% %             S_metaMin = S_metaMin + 2*pi^2/self.kappa*(Rho_metaMin(end)^2*dRho_metaMin(end));
% %         end
% %         
% %         % Add volume term for bubble interior
% %         volume = 8*pi^2/3/w^4*sin(R(1)*w/2)^4*(cos(R(1)*w)+2);
% %         S_metaMin = S_metaMin + volume*(2*self.V(self.phi_metaMin)) - ...
% %             16*pi^2*sin(R(1)*w/2)^4*(cos(R(1)*w)+2)/3/self.kappa/w^2;
% %         
% %         %% Compute tunneling rate exponent
% %         
% %         B = S_bubble-S_metaMin;
% %         
% %     end
% %     
% %     function B = find_exponent(self,R,Y)
% %         % Calculate the Euclidean action for the instanton:
% %         % 
% %         % .. math::
% %         %   S = \int [(d\phi/dr)^2 + V(\phi)] r^\alpha dr d\Omega_\alpha
% %         
% %         % f[r_?NumberQ] := - r^(rDim-1) profile[r] . (gradient[profile[r]]-gradient[falseVacuum]);
% %         % Pi^(rDim/2)/Gamma[1+rDim/2] NIntegrate[f[r], Join[{r, 0}, ri, {Infinity}],
% %         
% %         Y = mat2cell(Y,ones(1,size(Y,1)),size(Y,2));
% %         
% %         %% Get profile for false vacuum background
% %         
% %         Phi_metaMin  = self.phi_metaMin*ones(size(R));
% %         dPhi_metaMin = zeros(size(R));
% %         
% %         % Get Rho(R) in the case of uniform false vacuum
% %         w = abs(self.kappa*self.V(self.phi_metaMin))^0.5;
% %         if w == 0 || self.no_gravity
% %             Rho_metaMin   = R;
% %             dRho_metaMin  = ones(size(R));
% %             d2Rho_metaMin = zeros(size(R));
% %         elseif self.V(self.phi_metaMin) > 0
% %             Rho_metaMin   = sin(R*w)/w;
% %             dRho_metaMin  = cos(R*w);
% %             d2Rho_metaMin = -sin(R*w)*w;
% %         else
% %             Rho_metaMin   = sinh(R*w)/w;
% %             dRho_metaMin  = cosh(R*w);
% %             d2Rho_metaMin = sinh(R*w)*w;
% %         end
% %         
% %         Y_metaMin = {Phi_metaMin,dPhi_metaMin,Rho_metaMin,dRho_metaMin,d2Rho_metaMin};
% %         
% %         %% Compute tunneling rate exponent
% %         
% %         S_bubble  = self.find_action_4D(R,Y);
% %         S_metaMin = self.find_action_4D(R,Y_metaMin);
% %         
% %         B = S_bubble-S_metaMin;
% %         
% %     end
    


% %     function [r,y,convergence_type] = integrate_profile_old(self,r0,y0,dr0,epsabs,drmin,rmax)
% %         % Integrate the bubble wall equation
% %         % 
% %         % This works the same basic way as when there is no gravity, except
% %         % now the overshoot/undershoot conditions are a little bit more
% %         % complicated. The solution will overshoot if at any point `rho(r)
% %         % = 0` for `r > 0` (indicating that the bubble has wrapped around
% %         % to the anti-podal point of de Sitter space), unless `dphi/dr` is
% %         % also zero.
% %         
% %         dr = dr0;
% %         
% %         % dY is the ODE that we use
% %         dY = @(y,r) self.equation_of_motion(y,r);
% %         dydr0 = dY(y0,r0);
% %         
% %         ysign = sign(y0(1)-self.phi_metaMin); % positive means we're heading down, negative means heading up.
% %         rmax = rmax + r0;
% %         
% %         while true
% %             
% %             % Perform integration step
% %             options = odeset(...
% %                 'Jacobian', @(t,y) self.eom_jacobian(t,y),...
% %                 'Events',   @(t,y) self.myEvent(t,y,epsabs,rmax,ysign) );
% %             [r1,y1,re,ye,ie] = ode23s(@(r,y) dY(y,r).',[r0,r0+2],y0.',options);
% %             
% % %             if any(diff(r1) <= drmin/10)
% % %             disp(r1(end))
% %             if r1 == r0
% %                 drnext = 0.2*dr;
% % %                 ur = r1([true;diff(r1)/drmin > 10]);
% % %                 drnext = ur(end-1)-ur(1);
% % %                 drnext = min(y1(end,1)./y1(end,2));
% %                 dr = drnext;
% %                 continue
% %             else
% %                 drnext = 1.5*dr;
% % %                 drnext = min(y1(end,1)./y1(end,2));
% %             end
% %             
% %             r1 = r1(end,:);
% %             y1 = y1(end,:);
% %             
% %             dydr1 = dY(y1,r1);
% %             
% %             phi1 = y1(1); dphi1 = y1(2);
% %             rho1 = y1(3); drho1 = y1(4);
% %             
% %             % Check for errors
% %             if r1 > rmax
% %                 error('IntegrationError: r > rmax');
% %             elseif dr < drmin
% %                 % We need to check to see if we're at rho ~ 0. 
% %                 % If we are, then the equations are singular.
% %                 % But, it's a perfectly valid solution if dphi goes to 0 faster.
% %                 if drho1 < 0 && rho1 < 25*epsabs(3)
% %                     % linearly extrapolate to where rho = 0:
% %                     dr = -rho1/drho1;
% %                     r = r1 + dr;
% %                     y = y1 + dydr1*dr;
% %                     if abs(y(2)) < 3*epsabs(2)
% %                         % phi derivative also goes to zero
% %                         convergence_type = 'converged';
% %                     else
% %                         convergence_type = 'overshoot';
% %                     end
% %                     return
% %                 else
% %                     error('IntegrationError: dr < drmin'); % RNG
% %                 end
% %             end
% %             
% %             % Check for completion
% %             if abs(phi1-self.phi_metaMin) < 3*epsabs(1) && abs(dphi1) < 3*epsabs(2)
% %                 % Got close enough to meta-minimum with small enough dphi
% %                 
% %                 r = r1;
% %                 y = y1;
% %                 convergence_type = 'converged';
% %                 break
% %                 
% %             elseif dphi1*ysign > 0
% %                 % Didn't make it over the potential barrier, so the rate of
% %                 % change of phi has changed sign
% %                 
% %                 % Interpolate to where dphi(r) = 0
% %                 f = cubicInterpFunction(y0,dr*dydr0,y1,dr*dydr1);
% %                 
% %                 if feval(@(a) a(2), f(0).*f(1)) <= 0
% %                     x = fzero(@(x) feval(@(a) a(2),f(x)),[0 1]);
% %                     r = r0 + dr*x;
% %                     y = f(x);
% %                 else
% %                     r = r0;
% %                     y = nan(1,4);
% %                 end
% %                 
% %                 convergence_type = 'undershoot';
% %                 break
% %                 
% %             elseif (phi1-self.phi_metaMin)*ysign < 0
% %                 % Passed meta-minimum without stopping
% %                 
% %                 % Interpolate to where phi(r) = phi_metaMin
% %                 f = cubicInterpFunction(y0,dr*dydr0,y1,dr*dydr1);
% % 
% %                 if feval(@(a) a(1), (f(0)-self.phi_metaMin).*(f(1)-self.phi_metaMin)) < 1
% %                     x = fzero(@(x) feval(@(a) a(1),f(x))-self.phi_metaMin,[0,1]);
% %                     r = r0 + dr*x;
% %                     y = f(x);
% %                 else
% %                     r = r0;
% %                     y = nan(1,4);
% %                 end
% %                 
% %                 convergence_type = 'overshoot';
% %                 break
% %                 
% %             end
% %             
% %             % Advance the integration variables
% %             [r0,y0,dydr0] = deal(r1,y1,dydr1);
% %             dr = drnext;
% %             
% %         end
% %         
% %         % Check convergence for a second time. 
% %         % The extrapolation in overshoot/undershoot might have gotten us within
% %         % the acceptable error.
% %         if all(abs(y(1:2) - [self.phi_metaMin,0]) < 3*epsabs(1:2))
% %             convergence_type = 'converged';
% %         end
% %         
% %     end
    

% %     function [R,Y,Rerr] = integrate_and_save_profile(self,R,y0,dr,drmin)
% %         % Integrate the bubble profile, saving the output in an array.
% %         
% %         N = length(R); % Number of points
% %         
% %         r0 = R(1); % Initial radius
% %         
% %         % Field values [phi(R) dphi(R)]
% %         Y = zeros(N,length(y0));
% %         Y(1,:) = y0;
% %         
% %         % dY is the ODE that we use
% %         dY = @(y,r) self.equation_of_motion(y,r);
% %         dydr0 = dY(y0,r0);
% %         
% %         Rerr = NaN;
% %         
% %         i = 2;
% %         while i <= N
% %             
% %             % Perform integration step
% %             [r1,y1] = ode15s(@(r,y) dY(y,r).',[r0,r0+0.5*dr,r0+dr],y0.');
% %             
% %             if dr >= drmin
% %                 
% %                 if r1 == r0
% %                     drnext = 0.2*dr;
% %                     dr = drnext;
% %                     continue
% %                 else
% %                     drnext = 1.2*dr;
% %                 end
% %                 
% %                 r1 = r1(end,:);
% %                 y1 = y1(end,:);
% %                 
% %             else
% %                 
% %                 y1 = y0 + (y1(end,:)-y0)*drmin/dr;
% %                 dr = drmin;
% %                 drnext = drmin;
% %                 r1 = r1(end,:);
% %                 if isnan(Rerr)
% %                     Rerr = r1;
% %                 end
% %                 
% %             end
% %             
% %             dydr1 = dY(y1,r1);
% %             
% %             % Fill the arrays, if necessary
% %             % If this integration step took us beyond the next radius to
% %             % process, then fill in that value by interpolating between the
% %             % previous result and the current result.
% %             if r0 < R(i) && R(i) <= r1
% %                 f = cubicInterpFunction(y0, dr*dydr0, y1, dr*dydr1);
% %                 while (i <= N && r0 < R(i) && R(i) <= r1)
% %                     x = (R(i)-r0)/dr;
% %                     Y(i,:) = f(x);
% %                     i = i + 1;
% %                 end
% %             end
% %             
% %             % Advance the integration variables
% %             [r0,y0,dydr0] = deal(r1,y1,dydr1);
% %             dr = drnext;
% %             
% %         end
% %         
% %         % Add d2Rho to output (for computing action)
% %         dYdR  = dY(Y,R);
% %         d2Rho = dYdR(:,end);
% %         Y = [Y d2Rho].';
% %         
% %     end
    

% %         if self.useThinWallApproximation
% %             
% %             %% Get bubble wall tension
% %             
% %             [Phi,~,~,~,~] = deal(Y{:});
% %             
% %             % Bubble radius (thin wall approximation?)
% %             R0 = R(1);
% %             
% %             % Calculate surface tension of bubble wall
% %             % M = \int_{phi_0}^{phi_r} 2 \sqrt{V(phi) - V(phi_fv)} d\phi
% %             x = linspace(min(Phi),max(Phi),100);
% %             x = linspace(self.V(self.phi_metaMin),self.V(self.phi_absMin),100);
% %             sig = trapz(x,sqrt(2*(self.V(x)-self.V(self.phi_metaMin))));
% %             
% %             H = abs(self.kappa/3*self.V(self.phi_metaMin))^-0.5;
% %             
% %             %% Use thin wall approximation to compute bubble action
% %             
% %             DV = abs(self.V(self.phi_metaMin) - self.V(self.phi_absMin));
% %             R0 = (3*sig)/sqrt(9*H^2*sig^2 + DV^2);
% %             
% %             x = 6*pi*sig^2/DV; 
% %             y = DV/(self.V(self.phi_metaMin)+self.V(Phi(1)));
% %             z = sqrt(1 + 2*x*y + x^2);
% %             
% %             S = 27*pi^2*sig*4/DV^3*(...
% %                 ( (1 + x*y) - z ) / ...
% %                 ( x^2*(y^2 - 1)*z ) );
% %             
% % %             S2 = 2*pi^2/3*H^(-4)*(-2*DV + ...
% % %                 (9*H^2*sig^2 + 2*DV^2)/sqrt(9*H^2*sig^2 + DV^2) );
% %             
% %             %% Compute pre-factor
% %             
% %             % Calculate pre-factor sans Hubble factor exp(3Ht)
% %             % exp(zeta_R'(-2)) = exp(-zeta_R(3)/(2*pi)^2) = 0.97001
% %             lambda = -(S + 3/8/self.V(self.phi_metaMin)) + log(4*0.97001*R0^2*sig^2);
% %             
% %         else
            
% %     function [R_int,Y_int] = make_interior_points(self,R,delta_phi0,max_interior_pts,rmin)
% %         % Figure out how this works
% %         %
% %         % Usage:
% %         %
% %         % [R_int,Y_int] = make_interior_points(self,R,delta_phi0,max_interior_pts,rmin);
% %         % 
% %         % % Add internal points to profile
% %         % R = horzcat(fliplr(R_int),R);
% %         % Y = horzcat(fliplr(Y_int),Y);
% %         
% %         if nargin < 4 || isempty(max_interior_pts)
% %             max_interior_pts = floor(sqrt(length(R)));
% %         elseif max_interior_pts == 0
% %             R_int = zeros(1,0);
% %             Y_int = zeros(5,0);
% %             return
% %         end
% %         
% %         % Select grid of radii
% %         dx0 = R(2)-R(1);
% %         if R(1) / dx0 <= max_interior_pts
% %             n = ceil(R(1)/dx0);
% %             R_int = linspace(0,R(1),n+1);
% %             R_int = R_int(end:-1:1);
% %         else
% %             n = max_interior_pts;
% %             % R(1) = dx0 * (n + a*n*(n+1)/2)
% %             a = (R(1)/dx0 - n) * 2/(n*(n+1));
% %             N = n:-1:1;
% %             R_int = R(1) - dx0*(N + 0.5*a*N.*(N+1));
% %             R_int(1) = 0.0; % enforce this exactly
% %         end
% %         
% %         R_int = max(R_int,rmin);
% %         
% %         % Initial conditions in vicinity of bubble interior
% %         Phi_int     = zeros(size(R_int));
% %         Phi_int(1)  = self.phi_absMin + delta_phi0;
% %         dPhi_int    = zeros(size(R_int));
% %         dPhi_int(1) = 0;
% %         
% %         % Potential in vicinity of bubble interior
% %         V0   = self.V(Phi_int(1));
% %         dV0  = self.dV(Phi_int(1));
% %         d2V0 = self.d2V(Phi_int(1));
% %         
% %         w = abs(self.kappa/3*V0)^0.5;
% %         if w == 0
% %             % No curvature
% %             R_int2      = R_int;
% %             Rho_int     = R_int;
% %             dRho_int    = ones(size(R_int));
% %             d2Rho_int   = zeros(size(R_int));
% %         elseif V0 > 0
% %             % Positive curvature
% %             R_int2      = tan(w*R_int)/w;
% %             Rho_int     = sin(w*R_int)/w;
% %             dRho_int    = cos(w*R_int);
% %             d2Rho_int   = -w*sin(w*R_int);
% %         else
% %             % Negative curvature
% %             R_int2      = tanh(w*R_int)/w;
% %             Rho_int     = sinh(w*R_int)/w;
% %             dRho_int    = cosh(w*R_int);
% %             d2Rho_int   = w*sinh(w*R_int);
% %         end
% %         
% %         % Assume delta_phi0 is small, so the field in the bubble interior
% %         % is close to the absolute minimum, and thus the potential can be
% %         % approximated as quadratic, for which there is an exact solution
% %         for i = 1:length(R_int)
% %             [Phi_int(i),dPhi_int(i)] = self.exact_solution( ...
% %                 R_int2(i),Phi_int(1),dV0,d2V0);
% %         end
% %         
% %         Y_int = vertcat(Phi_int,dPhi_int,Rho_int,dRho_int,d2Rho_int);
% %         
% %     end
    
% %                         error('FalseVacuumInstanton:HawkingMoss',...
% %                             'Found min of abs(dphi/dr); no tunneling solution available.');