% jmd
% 8.14.12

farmConfig;

%% non-isozyme synth lethals

fid = fopen(strcat(mainpath, pathsep, 'farm_data', pathsep, 'allGenesGrowth.txt'));
C = textscan(fid, '%s %u %f', 'HeaderLines',1);
fclose(fid);

% haveEffect & grow
maybeSL = C{2}==1 & C{3}>grThresh;
nonessGenes=C{1}(maybeSL);
pairs = combnk(nonessGenes, 2);

% write out potential pairs to test
writeCell(pairs, strcat(mainpath, pathsep, 'farm_data', pathsep, 'potential_SLs.txt'), '%s %s', {'gene1' 'gene2'}, ' ');

npairs=size(pairs,1);
grVec=ones(npairs,1);

h = waitbar(0, 'Double gene deletion analysis in progress ...');
for i=1:npairs
    modelDelTmp = deleteModelGenes(limedNC, pairs(i,:));
    fbaSolTmp=optimizeCbModel(modelDelTmp, 'max');
    % this doesn't ever seem to return NaN, i think stat~=1 -> f=0
    grVec(i)=fbaSolTmp.f;
    waitbar(i/npairs,h);
    if mod(i,1000)==0
        i
    end
end
close(h);

writeCell(pairs(grVec<grThresh,:), strcat(mainpath, pathsep, 'farm_data', pathsep, 'nonIsozymeSynthLethals.txt'), '%s %s', {'gene1' 'gene2'}, ' ');

%% isozyme SLs

% -- get rxns catalyzed by 2 isozymes --
% get rxns catalyzed by 2 genes, could be can be complex or isozymes
% rxnGeneMat binary sparse, rxns by genes
twoGeneRxns=sum(limedNC.rxnGeneMat,2)==2;
% for isozymes, look for 'or' in grRules
orRxnsCell=regexp(limedNC.grRules, 'or');
orRxns=~cellfun(@isempty, orRxnsCell);
% remove spontaneous rxns
spontRxnsCell=regexp(limedNC.grRules, 's0001');
nonSpontRxns=cellfun(@isempty, spontRxnsCell);
% cross 2 gene rxns with 'or' rxns & remove spont rxns for 2-gene isozymes
isoRxns=limedNC.rxns(twoGeneRxns & orRxns & nonSpontRxns);

% -- check essentiality of isoRxns --
gr=zeros(length(isoRxns),1);
essIsoGenePairs={};
ub0=limedNC.ub;
for i=1:length(isoRxns)
    rxnTmpInd=strmatch(isoRxns(i), limedNC.rxns);
    limedNC.ub(rxnTmpInd)=0;
    op=optimizeCbModel(limedNC);
    gr(i,1)=op.f;
    
    % get gene pair
    if op.f<grThresh
        essIsoGenePairs=[essIsoGenePairs; limedNC.genes( find(limedNC.rxnGeneMat(rxnTmpInd,:)) )' ];
    end
    
    % re-initialize ub
    limedNC.ub=ub0;
end

% --write--
% this writes the isozyme SLs, albeit with two duplicates, which we removed
% in post-processing
writeCell(essIsoGenePairs, strcat(mainpath, pathsep, 'farm_data', pathsep, 'isozymeSLs.txt'), '%s %s', {'gene1' 'gene2'}, ' ');
