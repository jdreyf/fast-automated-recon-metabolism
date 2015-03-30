% jmd
% 5.9.13
% farmConfig.m

% these are the paths i use on my windows box and unix server, respectively
% set the path to the parent directory of farm & farm_data
path_windows='Z:\reconstruct\farm_gcode';
path_nonwin='/msc/neurospora/FBA';

%% set path based on windows vs. unix
os = getenv('OS');
if strcmp(os, 'Windows_NT')
    mainpath = path_windows;
    % windows path separator
    pathsep='\';
else
    mainpath = path_nonwin;
    % on unix server, i don't have a waitbar
    warning('off', 'MATLAB:waitbar:DeprecatedBehavior');
    pathsep='/';
end

cd(mainpath);

%% add paths
addpath(strcat(mainpath, pathsep, 'farm'));
addpath(strcat(mainpath, pathsep, 'farm', pathsep, 'utils'));
addpath(strcat(mainpath, pathsep, 'farm', pathsep, 'MD-FBA'));

%% set params
format long;
changeCobraSolverParams('LP','printLevel', 0);
changeCobraSolver('gurobi5', 'all');
% made lower bound 0.02 on nov 14, 2012. didn't want it higher, since then our model would be predicted not to grow on acetate.
grThresh=0.02;

suppUb=3;

%% read nc model
% cobra expects one letter compartment abbreviations: adding compList doesn't seem to help.
nc0 = readCbModel(strcat(mainpath, pathsep, 'farm_data', pathsep, 'nc10.xml'));
nc=compliantSbmlModel2ascii(nc0);

% need model.match for removeRxns() cobra function
nc.match=zeros(1, length(nc.rxns));

%% read mat files for train/test
cd(strcat(mainpath, pathsep, 'farm_data'));
load -mat nuts.mat
load -mat testEss.mat
load -mat testNoness.mat
load -mat testSupp.mat
load -mat trainEss.mat
load -mat trainNoness.mat
load -mat trainSupp.mat

%% get rxns that limed-fba shouldn't dilute
cd(strcat(mainpath, pathsep, 'farm'));
% biomass mets have 'Composition' in name, eg 'CellWallComposition'
compos=regexpi(nc.rxns,'Composition');
isBiomass=cellfun(@(x) ~isempty(x), compos);
nonDiluteRxns = [findNondiluteRxns(nc, false)' nc.rxns(isBiomass)'];

%% limed
% in past, chose some common mets to not dilute
%limedNC=limedModel(nc, 10^-4, nonDiluteRxns, {'PROTON', 'WATER', 'OH', 'OXYGEN-MOLECULE', 'CARBON-DIOXIDE'});
limedNC=limedModel(nc, 10^-4, nonDiluteRxns);
