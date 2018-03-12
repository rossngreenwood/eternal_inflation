classdef SRStatus < uint8
    
    enumeration
        
        ok  (0)         % Slow roll approximation is valid
        rho (1)         % Potential energy is negative; big crunch eminant
        eps (2)         % Slow roll condition \epsilon^2 << 1 is violated
        eta (3)         % Slow roll condition \abs(\eta) << 1 is violated
        localmin (4)    % Passed a local minimum of the potential
        neglambda (5)   % 
        rho_neglambda (6)         % Potential energy is negative; big crunch eminant
        eps_neglambda (7)         % Slow roll condition \epsilon^2 << 1 is violated
        eta_neglambda (8)         % Slow roll condition \abs(\eta) << 1 is violated
        localmin_neglambda (9)    % Passed a local minimum of the potential
        null (255)
        
    end
    
    methods
        
        function [tf] = isok(status)
            tf = (status == SRStatus.ok);
        end
        
    end
    
    methods (Static)
        
        function [status] = compute_status(V,Vp,Vpp,phi,Mpl)
            if     V(phi) < 0,                status = SRStatus.rho;
            elseif (Vp(phi)/V(phi)).^2/(16*pi/Mpl^2) > 1, status = SRStatus.eps;
            elseif abs(Vpp(phi)/V(phi))/(8*pi/Mpl^2) > 1, status = SRStatus.eta;
            else                              status = SRStatus.ok;     end
        end
        
    end
    
end