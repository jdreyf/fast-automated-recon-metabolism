% jmd
% 12.6.12

%% config
farmConfig;

%% compare methods on test + training set

% NaN growth implies FBA was infeasible, b/c couldn't meet ATP requirement
% so NaN growth corresponds to no-growth

% limed-FBA
predVsObs(limedNC, 'FBA', grThresh, trainEss, trainNoness, testEss, testNoness, trainSupp, testSupp, nuts, suppUb);

% fba
predVsObs(nc, 'FBA', grThresh, trainEss, trainNoness, testEss, testNoness, trainSupp, testSupp, nuts, suppUb);

% -- md-fba --
% requires tomlab cplex
% requires too much computer time to do everything together, so parallelize
% by running each of the 4 sets below separately 

model=nc; method='MDFBA'; verbose=true;

% train essentials (radford)
grTrEss=delGenesInMedia(model, trainEss, method, 2);
printAcc(grTrEss,grThresh,'train ess',0);

% train non-essentials (broad collection)
grTrNe = singleGeneDeletion2(model, method, trainNoness, verbose);
printAcc(grTrNe,grThresh,'train non-ess',1);

% test essentials (radford intersect borkovich)
grTeEss = singleGeneDeletion2(model, method, testEss(:,2), verbose);
printAcc(grTeEss,grThresh,'test ess',0);

% test non-essentials (borkovich po1 collection)
grTeNe = singleGeneDeletion2(model, method, testNoness, verbose);
printAcc(grTeNe,grThresh,'test non-ess',1);

%% comprehensive rescue prediction
% which supplements rescue which auxotroph mutants?
[suppAuxNames, rescueMat, supps, genes]=suppsByAuxotrophs(limedNC, trainSupp, testSupp, grThresh, suppUb);

%% single gene ko
[growth,grRateKO,grRateWT,hasEffect] = singleGeneDeletion2(limedNC, 'FBA');
limedEss=limedNC.genes(isnan(growth)|growth<grThresh);
% remove spontaneous rxn 'gene' s0001
limedEss=setdiff(limedEss, 's0001');
% write multiple fields from singleGeneDeletion2
fileID = fopen(strcat(mainpath, pathsep, 'farm_data', pathsep, 'allGenesGrowth.txt'), 'w');
fprintf(fileID,'%s %s %s\n', 'gene', 'hasEffect', 'growth');
for i=1:length(limedNC.genes)
    fprintf(fileID,'%s %u %f\n', limedNC.genes{i}, hasEffect(i), growth(i));
end
fclose(fileID);

%% comprehensive synthetic lethals
% see allSynthLethals.m