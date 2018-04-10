%% Estructuracion del programa

% 1 - Carga del phantom/paciente
% 2 - Carga de parametros alpha y beta
% 3 - Definicion de restricciones "distintas"
% 4 - Introduccion de los datos basicos
% 5 - Cálculos

%% 1 - Carga del phantom/paciente

close all
clear
if ~ismac
    opengl software
end

load HEAD_AND_NECK.mat; phantomtype = 'Head and Neck';

%% 2 - Carga de parametros alpha y beta

cst = prueba_abLoader (cst, phantomtype);

%% 3 - Definicion de restricciones 

%       .type                   || NaN parameters

% 'square underdosing'          ||.EUD - .volume
% 'square overdosing'           ||.EUD - .volume
% 'square deviation'            ||.EUD - .volume
% 'mean'                        ||.dose - .EUD - .volume
% 'EUD'                         ||.dose - .volume
% 'min dose constraint'         ||.penalty - .EUD - .volume
% 'max dose constraint'         ||.penalty - .EUD - .volume
% 'min mean dose constraint'    ||.EUD - .volume
% 'max mean dose constraint'    ||.EUD - .volume
% 'min EUD constraint'          ||.penalty - .dose - .volume
% 'max EUD constraint'          ||.penalty - .dose - .volume
% 'min DVH constraint'          ||.penalty - .EUD 
% 'max DVH constraint'          ||.penalty - .EUD
% 'min DVH objective'           ||.EUD
% 'max DVH objective'           ||.EUD
    
% Parotid_LT
cst{13,6}.type = 'square overdosing';
cst{13,6}.dose = 25;
cst{13,6}.EUD = NaN;
cst{13,6}.penalty = 100;
cst{13,6}.volume = NaN;
cst{13,6}.robustness = 'none';

% Parotid_RT
cst{14,6}.type = 'square overdosing';
cst{14,6}.dose = 25;
cst{14,6}.EUD = NaN;
cst{14,6}.penalty = 100;
cst{14,6}.volume = NaN;
cst{14,6}.robustness = 'none';

% PTV63
cst{15,6}.type = 'square deviation';
cst{15,6}.dose = 63;
cst{15,6}.EUD = NaN;
cst{15,6}.penalty = 1000;
cst{15,6}.volume = NaN;
cst{15,6}.robustness = 'none';

% PTV70
cst{16,6}.type = 'square deviation';
cst{16,6}.dose = 63;
cst{16,6}.EUD = NaN;
cst{16,6}.penalty = 1000;
cst{16,6}.volume = NaN;
cst{16,6}.robustness = 'none';

% SKIN
cst{17,6}.type = 'square overdosing';
cst{17,6}.dose = 30;
cst{17,6}.EUD = NaN;
cst{17,6}.penalty = 800;
cst{17,6}.volume = NaN;
cst{17,6}.robustness = 'none';


%% 4 - Introduccion de los datos basicos

% meta information for treatment plan 
pln.numOfFractions = 30;
pln.radiationMode = 'protons';           % either photons / protons / helium / carbon
pln.machine = 'Generic';

% beam geometry settings
pln.propStf.bixelWidth     = 5; % [mm] / also corresponds to lateral spot spacing for particles
pln.propStf.gantryAngles   = [0]; % [?];
pln.propStf.couchAngles    = [0]; % [?];
pln.propStf.numOfBeams     = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter      = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% optimization settings
pln.propOpt.runDAO         = false;   % 1/true: run DAO, 0/false: don't / will be ignored for particles
pln.propOpt.runSequencing  = false;   % 1/true: run sequencing, 0/false: don't / will be ignored for particles and also triggered by runDAO below


%% 5 - Cálculos

[~, ResultConstRBE, ResultRBEMCN, ResultRBEUCM, DoseStatistics, NTCP] = prueba_NTCP(cst, pln, ct, phantomtype);
