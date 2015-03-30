% jmd
% 6.25.12
% delGenesInMedia.m

function [growth]=delGenesInMedia(model, condCell, method, geneCol, supps, suppUb)
if nargin<3
    method='FBA';
end
if nargin<4
    geneCol=2;
end
if nargin<5
    supps={};
end
if nargin<6
    suppUb=10;
end

%index of column of cpds to remove from medium to create supplemental media
rmCpdCol=geneCol+1;
%index of column of cpds to add to medium to create supplemental media
addCpdCol=rmCpdCol+1;

%initialize growth vector to NaN
growth=nan(size(condCell,1),1);

for i=1:size(condCell,1)
    model2=model;
    if ~isempty(condCell{i, rmCpdCol})
        rmRxnName=strcat(condCell{i, rmCpdCol}, '-TRANS-RXN-L2R');
        model2=removeRxns(model2, rmRxnName, true, false);
    end
    
    % add rxn from addCpdCol
    if ~isempty(condCell{i, addCpdCol}) 
        addMetName=strcat(condCell{i, addCpdCol}, '[CCO-EXTRACELLULAR]');
        addRxnName=strcat(condCell{i, addCpdCol}, '-TRANS-RXN-L2R');
        model2=addReaction(model2, addRxnName, {addMetName}, 1, false, 0, suppUb); 
    end
    
    if ~isempty(supps)
        % suppExpr may be "MET and THR"
        suppTmpExpr=splitString(supps{i}, ' ');
        suppTmp=setdiff(suppTmpExpr, 'and');
        % loop thru supps in supps{i}
        for suppTmpInd=1:length(suppTmp)
            addMetName=strcat(suppTmp{suppTmpInd}, '[CCO-EXTRACELLULAR]');
            addRxnName=strcat(suppTmp{suppTmpInd}, '-TRANS-RXN-L2R');
            % if rxn of same name already exists, don't add again
            if ~ismember(addRxnName, model2.rxns)
                model2=addReaction(model2, addRxnName, {addMetName}, 1, false, 0, suppUb);
            end
        end
    end
    
    geneKOs=condCell{i,geneCol};
    if ismember(geneKOs, model2.genes)
        model2 = deleteModelGenes(model2, geneKOs);
    end
    
    if strcmpi(method,'FBA')
        sol=optimizeCbModel(model2);
        growth(i,1)=sol.f;
    else
        sol=mdFBA(model2,0);
        growth(i,1)=sol.result_opt;
    end
end