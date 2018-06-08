%% 
close all
clear
if ~ismac
    opengl software
end
load PROSTATE.mat; phantomtype = 'Prostate';


%% 

% % PTV_68
% cst{6,5}.Priority = 1;
% cst{6,6}.type = 'max dose constraint';
% cst{6,6}.dose = 68;
% cst{6,6}.penalty = NaN;
% cst{6,6}.EUD = NaN;
% cst{6,6}.volume = NaN;
% cst{6,6}.robustness = 'none';

%%
% PTV_68
cst{8,5}.Priority = 3;
cst{8,6}.type = 'max dose constraint';
cst{8,6}.dose = 50;
cst{8,6}.penalty = NaN;
cst{8,6}.EUD = NaN;
cst{8,6}.volume = NaN;
cst{8,6}.robustness = 'none';

%%
% Rectum
% cst{1,5}.Priority = 2;
% cst{1,6}.type = 'max DVH constraint';
% cst{1,6}.dose = 50;
% cst{1,6}.penalty = NaN;
% cst{1,6}.EUD = NaN;
% cst{1,6}.volume = 1;
% cst{1,6}.robustness = 'none';
%%
% PTV_56
% cst{1,5}.Priority = 2;
% cst{1,6}.type = 'max DVH constraint';
% cst{1,6}.dose = 56;
% cst{1,6}.penalty = NaN;
% cst{1,6}.EUD = NaN;
% cst{1,6}.volume = 1;
% cst{1,6}.robustness = 'none';


%% 

% meta information for treatment plan
pln.numOfFractions = 39;
pln.radiationMode = 'protons';           % either photons / protons / helium / carbon
pln.machine = 'Generic';

% beam geometry settings
pln.propStf.bixelWidth     = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles   = [90 270]; % [?];
pln.propStf.couchAngles    = [0 0]; % [?];
pln.propStf.numOfBeams     = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter      = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln.propOpt.runDAO         = false;   % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing  = false;   % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below


quantityOpt         = 'RBExD';
modelName           = 'MCN';             % none: for photons, protons, carbon                                    MCN: McNamara-variable RBE model for protons
% WED: Wedenberg-variable RBE model for protons

scenGenType = 'nomScen';            % scenario creation type 'nomScen' 'wcScen' 'impScen' 'rndScen'

% retrieve model parameters
pln.bioParam = matRad_bioModel(pln.radiationMode,quantityOpt, modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_multScen(ct,scenGenType);
stf = matRad_generateStf(ct,cst,pln);

     
%% dose calculation
if strcmp(pln.radiationMode,'photons')
    dij = matRad_calcPhotonDose(ct,stf,pln,cst);
elseif strcmp(pln.radiationMode,'protons') || strcmp(pln.radiationMode,'helium') || strcmp(pln.radiationMode,'carbon')
    dij = matRad_calcParticleDose(ct,stf,pln,cst);
end

resultGUI = matRad_fluenceOptimization(dij,cst,pln);
pause(1.5);
close all
clc


matRadGUI