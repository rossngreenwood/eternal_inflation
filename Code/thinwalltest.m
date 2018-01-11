
xguess           = [];
xtol             = 1e-4;
phitol           = 1e-4;
thinCutoff       = 1e-4;
rmin             = 1e-4;
rmax             = 1e+4;

mv = 1e0;
mh = 0.45205;
mh = 0.0001;
mh = 1;

% False and true vacua
phi_metaMin = -0.55*mh;
phi_absMin  = +7.45*mh;

% (3*(x/mh) - (x/mh).^3 + .2*(x/mh).^4 + 10)*mv
% ((x/mh).^2/2 - (x/mh).^3/2 + a*(x/mh).^4/8)*mv
a = 0.99;

% Initialize object
fvi = FalseVacuumInstanton(...
    'V',            @(x) (3*(x/mh) - (x/mh).^3 + .2*(x/mh).^4 + 10)*mv,...
    'dV',           @(x) (3 - 3*(x/mh).^2 + .8*(x/mh).^3)*mv/mh,...
    'd2V',          @(x) (- 6*(x/mh) + 2.4*(x/mh).^2)*mv/mh^2,...
    'M_Pl',         sqrt(8*pi),...
    'phi_metaMin',  phi_metaMin,...
    'phi_absMin',   phi_absMin,...
    'no_gravity',   false );

tic
% Get profile
[R,Y,useThinWall] = fvi.find_profile(xguess,xtol,phitol,thinCutoff,rmin,rmax);
toc

B = fvi.find_tunneling_suppression(R,Y)