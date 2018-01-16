% False and true vacua
phi_metaMin = +3;
phi_absMin  = -2;

xguess           = [];
xtol             = 1e-4;
phitol           = 1e-4;
thinCutoff       = 1e-4;
rmin             = 1e-4;
rmax             = 1e+4;

% Initialize object
fvi = FalseVacuumInstanton(...
    'V',            @(x) (x.^4 - 14*x.^2 + 24*x)*1e-0,...
    'dV',           @(x) (4*x.^3 - 28*x + 24)*1e-0,...
    'd2V',          @(x) (12*x.^2 - 28)*1e-0,...
    'M_Pl',         500,...
    'phi_metaMin',  phi_metaMin,...
    'phi_absMin',   phi_absMin,...
    'no_gravity',   true);

tic
% Get profile
[R,Y] = fvi.find_profile(xguess,xtol,phitol,thinCutoff,rmin,rmax);
toc

B = fvi.find_tunneling_suppression(R,Y)
% B = -log(lambda);

plot(R,Y(:,1))