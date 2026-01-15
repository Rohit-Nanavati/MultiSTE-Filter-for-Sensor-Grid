% clc;
close all;
clear;
set(0,'defaultAxesFontSize',18);
set(0,'defaultAxesFontName','Times New Roman');
set(0,'defaultTextFontName','Times New Roman');
set(0,'DefaultFigureColormap',flipud(hot));

addpath('Functions');
% rng('default')

plotsethandle = @mste_figplot_set;
evn_scaling = 1;
%% Simulated source parameters and Ground Truth Generation

% True source
ts.Q = [6 7]*evn_scaling;
ts.x = [10 40]*evn_scaling;
ts.y = [40 30]*evn_scaling;

ts.u = 4; % wind speed
ts.phi = deg2rad(-90);
ts.ci = 1.2;% Also s.D for Pasquil model
ts.cii = 5;% Also s.t or tau for Pasquil Model


% for particle filter
N_sources = 4; %numel(ts.Q);

D_threshold = 5e-4;
windmax = ts.u;
what = ts.u/windmax;
windstruct.ratio = what;
windstruct.heading = wrapToPi(ts.phi);

PDetect = 0.95;
PMissDetect = 1-PDetect;

% Create rectangular domain area
xlimits = [0 50]*evn_scaling;
ylimits = [0 50]*evn_scaling;
xmin = xlimits(1); xmax = xlimits(2);
ymin = ylimits(1); ymax = ylimits(2);

InpDom = [xlimits([1 2 2 1]);
    ylimits([1 1 2 2])]';

Domain = polyshape(InpDom);

% example data
xstart = 5; xend = 45;
ystart = 5; yend = 45;
NofColumns = 5;
gridresolution = [(xend-xstart) (yend-ystart)]/(NofColumns-1);
xpose_vec = xstart:gridresolution(1):xend;
ypose_vec = ystart:gridresolution(2):yend;

[x_matrix,y_matrix] = meshgrid(xpose_vec,ypose_vec);
pose_vec = [x_matrix(:) y_matrix(:)];

%% Initialize Particle Filter (PF)
ParticleFilterFunction = @PFupdateSimple;
TransitionSigma.x = 1;
TransitionSigma.y = 0.5;
TransitionSigma.Q = 0.3;


GMMixCompNo_max = ceil(N_sources);
Np = 25000; % No of particle for the PF
theta.x = zeros(Np, N_sources); % Each particle has x-coordinates for all sources
theta.y = zeros(Np, N_sources); % y-coordinates for all sources
theta.Q = zeros(Np, N_sources); % Release rates for all sources

Kmax = N_sources;            % keep your config var but treat as max
theta.mask = ones(Np, Kmax);    % start with all active, or set a prior pattern

for ns = 1:N_sources
    % Uniform prior for location
    theta.x(:,ns) = xmin + (xmax-xmin) * rand(Np,1);
    theta.y(:,ns) = ymin + (ymax-ymin) * rand(Np,1);

    % Gamma Distribution for the release rate
    a = ones(Np,1)*2;
    b = ones(Np,1)*5;
    theta.Q(:,ns) = gamrnd(a,b);

    % theta.Q(:,ns) = ts.Q(ns) + 4*randn(Np,1);
end

% ----- Equal prior over M: Uniform{1..Kmax} -----
Mrand = randi([1, Kmax], Np, 1);         % each particle gets an M_i
theta.mask = zeros(Np, Kmax);        % 0=inactive, 1=active
for i = 1:Np
    cols = randperm(Kmax, Mrand(i));     % pick which sources are active
    theta.mask(i, cols) = 1;
end

% Enforce inactivity: Q=0 where mask==0 (so inactive sources contribute nothing)
theta.Q = theta.Q .* theta.mask;

% Particle weights initialization with equal weights
Wp = ones(Np,1);
Wpnorm = Wp./sum(Wp);

%% Loop parameters
t_now = 0;
itermax = 600;
filter_counter_max = 35;
normed_std = zeros(2,itermax);
ste_esterror = zeros(1,itermax);
stePose_esterror = zeros(1,itermax);
telapsed_Density = zeros(1,itermax);
telapsed_control = zeros(1,itermax);
telapsed_PF = zeros(1,itermax);

iter = 1;
percentageNoise = 0.3;
% Initially sampling the environment for preliminary belief update
EstHistory(iter).pose = pose_vec;
NoS = size(pose_vec,1);
D_sim = Plume_model(ts,pose_vec,ts);
D_measure = D_sim + percentageNoise*D_sim.*randn(NoS,1);
D_measure(D_measure < D_threshold) = 0;
D_measure(rand(NoS,1)<=PMissDetect) = 0;
EstHistory(iter).measurement = D_measure;
fprintf('-Initial Location Sampled-\n')

%% Particle Filter Belief Update
ClusterTransitionFlag = 1;
fprintf("-Filter update no.: 1-\n");
tStart_PF = tic;
[theta, Wpnorm] = ParticleFilterFunction(theta,Wpnorm,D_measure,pose_vec,Np,D_threshold,ts,PDetect,TransitionSigma,Domain);
telapsed_PF(iter) = toc(tStart_PF);
filter_counter = 1;
[EstHistory(iter).confidentest,EstHistory(iter).estconfidence,~,~,~,EstHistory(iter).scatterplotpoints] = EstimateCompute(theta,Wpnorm);

EstStruct = EstHistory(iter).confidentest;
GOSPA = MultiSourceSTE_GOSPA(EstStruct.mean,ts);
normed_std(1,iter) = EstStruct.normedstd;
normed_std(2,iter) = EstStruct.normedPosestd;
ste_esterror(1,iter) = GOSPA.EstError;  
stePose_esterror(1,iter) = GOSPA.EstPoseError;   

tsamplingtime = 5;
tnow = 0;
tvec(iter) = tnow;

EstHistory(iter).particles = theta;
EstHistory(iter).weights = Wpnorm;
[f,at1,at2,plotstruct] = mste_figplot_initial(EstHistory,Domain,ts,1,tnow,N_sources);
OptimizationCost = zeros(NoS,itermax);

while iter<=itermax % normed_std(1,iter) > sigmaThreshold && iter<=itermax

    % Plotting function
    [f,at1,at2,plotstruct] = plotsethandle(EstHistory,theta,tnow,f,at1,at2,plotstruct);

    fprintf("--Sampling Environment--\n\n");

    D_sim = Plume_model(ts,pose_vec,ts);
    D_measure = D_sim + percentageNoise*D_sim.*randn(NoS,1);
    D_measure(D_measure < D_threshold) = 0;
    D_measure(rand(NoS,1)<=PMissDetect) = 0;
    EstHistory(iter).measurement = D_measure;

    % Particle Filter Belief Update
    fprintf("-Filter update no.: %d-\n",filter_counter+1);
    filter_counter = filter_counter +1;
    tStart_PF = tic;
    [theta, Wpnorm] = ParticleFilterFunction(theta,Wpnorm,D_measure,pose_vec,Np,D_threshold,ts,PDetect,TransitionSigma,Domain);
    telapsed_PF(filter_counter) = toc(tStart_PF);
    tnow = tnow + tsamplingtime;
    
    iter = iter + 1; % Increment iteration counter
    tvec(iter) = tnow;

    [EstHistory(iter).confidentest,EstHistory(iter).estconfidence,~,~,~,EstHistory(iter).scatterplotpoints] = EstimateCompute(theta,Wpnorm);

    % Evaluation Metrics
    EstStruct = EstHistory(iter).confidentest;
    GOSPA = MultiSourceSTE_GOSPA(EstStruct.mean,ts);
    normed_std(1,iter) = EstStruct.normedstd;
    normed_std(2,iter) = EstStruct.normedPosestd;
    ste_esterror(1,iter) = GOSPA.EstError;   
    stePose_esterror(1,iter) = GOSPA.EstPoseError;
   
    EstHistory(iter).particles = theta;
    EstHistory(iter).weights = Wpnorm;

end

% Final Plot Update
[f,at1,at2,plotstruct] = plotsethandle(EstHistory,theta,tnow,f,at1,at2,plotstruct);
% 
% figure(2);
% plot(tvec,stePose_esterror(1:iter)); hold on;
% plot(tvec,normed_std(2,1:iter));
% xlabel('Time'); 
% legend('Localisation GOSPA','Localisation Uncertainity');
% 
% figure(3);
% plot(tvec,ste_esterror(1:iter)); hold on;
% plot(tvec,normed_std(1,1:iter));
% yline(sigmaThreshold,'LineWidth',2);
% xlabel('Time'); 
% legend('GOSPA','Estimation Uncertainity');