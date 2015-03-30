% jmd
% 8.11.12
% loop thru all genes x all supps to see what grows

% init mat, loop thru genes, assign vector of supps, checkGrowth, report
function [suppAuxNames, rescueMat, supps, genes] = suppsByAuxotrophs(model, trainSupp, testSupp, grThresh, suppUb, method)
if nargin<6
    method='FBA';
end

%combine train+test gene-by-supps sets
suppCellAll=[trainSupp' testSupp']';
supps=unique(suppCellAll(:,5));
%get rid of genes (NCU's) incorrectly predicted as non-essential
predNeGenes=delGenesInMedia(model, suppCellAll, method);
ncuGrow=suppCellAll(~isnan(predNeGenes) & predNeGenes>grThresh, 2);
suppCellAll=suppCellAll(~ismember(suppCellAll(:,2), ncuGrow),:);
% subset to unique genes
% some genes associated w/ mult rescues
% each row a gene is on has same add+rm Cpd, so this is retained
[auxs,uGenesInd]=unique(suppCellAll(:,2));
suppCell=suppCellAll(uGenesInd,:);
% sort
[B,rowInd] = sort(suppCell(:,1));
suppCell=suppCell(rowInd,:);
% get dimensions
nsupp=length(supps);
naux=size(suppCell,1);
genes=suppCell(:,1);

% create rescue matrix
rescueMat=zeros(nsupp,naux);
for suppTmpInd=1:nsupp
    % make vector of supp
    suppTmpCell=cell(naux,1);
    suppTmpCell(1:naux,1)=supps(suppTmpInd);
    
    rescueMat(suppTmpInd,:)=delGenesInMedia(model, suppCell(:,1:4), method, 2, suppTmpCell, suppUb);
    fprintf('Finished %d out of %d supps\n', suppTmpInd, nsupp);
end

% create cell array of (gene,supp) name pairs
suppAuxNames={};
for i=1:nsupp
    for j=1:naux
        if rescueMat(i,j)>grThresh
            suppAuxNames=[{suppCell(j,1) supps(i)}', suppAuxNames];
        end
    end
end
suppAuxNames=suppAuxNames';