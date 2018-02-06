% False and true vacua
phi_metaMin = 2;
phi_absMin  = -3;

xguess           = [];
xtol             = 1e-4;
phitol           = 1e-4;
thinCutoff       = 1e-4;
rmin             = 1e-4;
rmax             = 1e+4;

factor = sqrt(10^randi([-4,4]));
factor = 1;
log10(factor^2)

% Initialize object
fvi = FalseVacuumInstanton(...
    'V',            @(x) (x.^4 - 14*x.^2 + 24*x)*factor^2,...
    'dV',           @(x) (4*x.^3 - 28*x + 24)*factor^2,...
    'd2V',          @(x) (12*x.^2 - 28)*factor^2,...
    'M_Pl',         sqrt(8*pi),...
    'phi_metaMin',  phi_metaMin,...
    'phi_absMin',   phi_absMin,...
    'no_gravity',   true);

tic
% Get profile
[R,Y] = fvi.find_profile(xguess,xtol,phitol,thinCutoff,rmin/sqrt(10),rmax/sqrt(10));
toc

% size(Y)
% Y(:,2) = Y(2)*factor;
% Y(:,4) = Y(4)*factor;
% Y(:,5) = Y(5)*factor^2;

B = fvi.find_tunneling_suppression(R,Y)*factor^2
% B = -log(lambda);

plot(R,Y(:,1))