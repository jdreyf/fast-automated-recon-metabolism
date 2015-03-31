% jmd
% 8.7.12

function []=predVsObs(model, method, grThresh, trainEss, trainNoness, testEss, testNoness, trainSupp, testSupp, nuts, suppUb)
if nargin<2
    method='FBA';
end
if nargin<3
    grThresh=10^(-6);
end
if ~ismember(method, {'FBA', 'MDFBA'})
    warning('fba:method', 'method not one of FBA or MDFBA');
end

if strcmp(method, 'MDFBA')
    verbose=true;
else
    verbose=false;
end

%% train essentials (Radford e-Compendium)
grTrEss=delGenesInMedia(model, trainEss, method, 2);
printAcc(grTrEss,grThresh,'train ess',0);

%% train non-essentials (broad collection)
% we get 1 wrong from 110, fba none wrong
grTrNe = singleGeneDeletion2(model, method, trainNoness, verbose);
printAcc(grTrNe,grThresh,'train non-ess',1);

%% test essentials (radford intersect borkovich)
grTeEss = singleGeneDeletion2(model, method, testEss(:,2), verbose);
printAcc(grTeEss,grThresh,'test ess',0);

%% test non-essentials (borkovich po1 collection)
grTeNe = singleGeneDeletion2(model, method, testNoness, verbose);
printAcc(grTeNe,grThresh,'test non-ess',1);

%% nuts
gr=delGenesInMedia(model, nuts, method, 2);
printAcc(gr,grThresh,'essential nuts',0);

gr=delGenesInMedia(model, nuts([1 4],:), method, 2, nuts([1 4], 5));
printAcc(gr,grThresh,'other media',1);

%% subset supps
ncuNoGr=[trainEss(isnan(grTrEss) | grTrEss<=grThresh, 2)' testEss(isnan(grTeEss) | grTeEss<=grThresh, 2)'];
trainSuppSS=trainSupp(ismember(trainSupp(:,2), ncuNoGr),:);
testSuppSS=testSupp(ismember(testSupp(:,2), ncuNoGr),:);

%% train supps
grTrSupps=delGenesInMedia(model, trainSuppSS, method, 2, trainSuppSS(:,5), suppUb);
printAcc(grTrSupps,grThresh,'train supp',1);

%% test supps
grTeSupps=delGenesInMedia(model, testSuppSS, method, 2, testSuppSS(:,5), suppUb);
printAcc(grTeSupps,grThresh,'test supp',1);