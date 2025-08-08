close all
clear all

path_to_MCX='/autofs/space/oxm_002/users/MCX19Workshop_DoNotModify/Software/linux64/MCXStudio'; %% PUT YOUR OWN PATH HERE

savepwd=pwd;
cd(path_to_MCX);
mcxsuite_addpath;
cd(savepwd);

rand_seed = randi([1 2^31-1], 1, 1);

%% Define optical parameters

%== Define simulated photon count (suggested less than 1e8 to avoid lengthy simulations)
cfg.nphoton=2e9;

%== Define simulated domain volume (3D array of labels or floating point numbers)
% define a 1cm radius sphere within a 6x6x6 cm box with a 0.5mm resolution
dimx=200;
dimy=200;
dimz=60; % in mm
[xi,yi,zi]=meshgrid(1:dimx,1:dimy,1:dimz);

cfg.vol=3*ones(size(xi)); % set to deep medium index
cfg.unitinmm=1; % define pixel size in terms of mm

scalp_th=5; idx_z_scalp=1:floor(scalp_th/cfg.unitinmm);
skull_th=7; idx_z_skull=(1+floor(scalp_th/cfg.unitinmm)):((scalp_th+skull_th)/cfg.unitinmm);

cfg.vol(:,:,idx_z_scalp)=1;
cfg.vol(:,:,idx_z_skull)=2;

%== Define optical properties for each tissue label
%         mua(1/mm) mus(1/mm)  g    n

cfg.prop=[0 0 1 1          % medium 0: the environment
    0.0164 0.74 0.01 1.37     % medium 1: scalp
    0.0115 0.81 0.01 1.37     % medium 1: skull
    0.017 1.16 0.01 1.37];   % medium 2: brain
	
	
%== Define source position (in grid-unit, not in mm unit!)
cfg.srcpos=[100,80,1];

%== Define source direction (a unitary vector)
cfg.srcdir=[0 0 1];

%== Define time-gate settings
cfg.tstart=0;   % starting time: 0
cfg.tend=1e-8;  % ending time: 10 ns
cfg.tstep=1e-8; % time step: 10 ns
cfg.seed = rand_seed;

%% Define detectors

%== Define detector position and radius
cfg.detpos=[100 105 1 1; 100 110 1 1.5;100 115 1 2];
% all options: 'dspmxvw'
% for details please type: help mcxlab

%== Define output structure
% d: detector id; s: scattering event count; p: partial path length
% x: exit position; v: exit direction
cfg.savedetflag = 'dspm';

numdet=size(cfg.detpos,1);
for det_idx=1:numdet,
    sdsep(det_idx)=norm(cfg.srcpos-cfg.detpos(det_idx,1:3));
end

%% Define GPU parameters

cfg.gpuid=1;         % use the first GPU
%cfg.gpuid='11';    % use two GPUs together
cfg.autopilot=1;     % assign the thread and block automatically
cfg.isreflect=1; % enable reflection at exterior boundary
cfg.debuglevel='P';

%% Preview the domain

mcxpreview(cfg);

%% Run simulation

[flux,detpos]=mcxlab(cfg);

%% generate g2s

mtau=logspace(-7,-1,500); % correlation lag times

disp_model='brownian';
lambda=850;
assumed_beta=0.5;

DV=[1e-6 1e-8 5e-6]; % BFi for scalp, skull, brain in units of mm^2/savedetflag

[mtau,g1]=generate_g1_mcxlab(cfg,detpos,mtau,disp_model,DV,lambda);
assumed_beta=0.5;

g2=1+assumed_beta*(g1.^2);


%% fit using semi-infinite analytical model

fit_options.lambda_dcs = 850*1e-6; % mm-1
fit_options.n=1.37;
fit_options.mu_a = 0.01; % mm-1
fit_options.mu_s = .8; % mm-1
fit_options.alpha = 1;
x0 = [0.5,2]; % beta, then Db times 1e9
lb = zeros(size(x0));
ub=[];

ft=1; lt=size(g2,2);   % could choose to fit less of the curve here

for detidx=1:numdet
    fit_options.rho=sdsep(detidx);
    test_x(detidx,:) = lsqcurvefit(@(x,taus)semi_infinite_g2(x,mtau(ft:lt),fit_options),x0,mtau(ft:lt),g2(detidx,ft:lt)',lb,ub);
end

fit_beta=test_x(:,1); 
fit_BFi=test_x(:,2)/1e9;






