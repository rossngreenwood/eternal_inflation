% Trial potential function
V = @(x) (x).^4 - 20*(x).^2 + 24*(x);

% False and true vacua
phi_metaMin = -2;
phi_absMin  = 3;

xguess           = 2;   
xtol             = 1e-4; 
phitol           = 1e-4; 
thinCutoff       = 1e-2; 
npoints          = 500;  
rmin             = 1e-5; 
rmax             = 1e4;  
max_interior_pts = 0;  

% profile on

% Initialize object
tic
fvi = FalseVacuumInstanton(...
    'V',@(x) (x.^4 - 14*x.^2 - 24*x)*1e-0,...
    'dV',@(x) (4*x.^3 - 28*x - 24)*1e-0,...
    'd2V',@(x) (12*x.^2 - 28)*1e-0,...
    'M_Pl',sqrt(8*pi).^3,...
    'phi_metaMin',phi_metaMin,...
    'phi_absMin',phi_absMin);

% Get profile
[R,Y,~,n_interior_pts] = fvi.find_profile(xguess,xtol,phitol,thinCutoff,npoints,rmin,rmax,max_interior_pts);

% Calculate Euclidean action
S = fvi.find_action(R,Y);
lambda = fvi.find_tunneling_rate(R,Y);

toc

plot(R,Y(1,:))

% profile viewer