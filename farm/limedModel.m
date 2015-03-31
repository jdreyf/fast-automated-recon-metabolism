% need to comment at top, make sure model is irreversible
% last 2 parameters are cells of names
function [model2] = limedModel(model, eps, nondiluteRxns, nonDiluteBaseMetNames)
%prepare S for linear metabolite dilution FBA (LIMED-FBA)

% need irrev model
if modelIsReversible(model)
    error('limedFBA:rev', 'an irreversible model is required');
end

if nargin < 2
    eps=10^-4;
end
if nargin < 3
    nondiluteRxns = {};
end
if nargin < 4
	nonDiluteBaseMetNames={};
end	

baseMetNames=arrayfun(@parseMetNames, model.mets);

%find indices of non-zero stoichiometric elements to potentially dilute
[nzRow, nzCol] = find(model.S~=0);
nzMat = [nzRow nzCol];

%get cpds to not dilute in any compartment
nondiluteMetsInd = find(ismember(baseMetNames, nonDiluteBaseMetNames));

%get rxns to not dilute
nondiluteRxnsInd = find(ismember(model.rxns, nondiluteRxns));

%set elements to dilute
%each row of nzMat corresponds to an element of model.S
nzMat = nzMat(~ismember(nzRow, nondiluteMetsInd) & ~ismember(nzCol, nondiluteRxnsInd),:);

%get sum of rxns diluted per metabolite
subsetS = model.S(:,setdiff(1:size(model.S,2), nondiluteRxnsInd));
numNonzeroRowSums = sum(subsetS~=0,2);

%metab conc go up to 10^-2 (bennet et al, nat chem bio, 2009) and we use ub = 10^3, so multiply by 10^-5
%this sometimes gives numerical error, though, so use 10^-4 and verify KOs using FBA: 
%limed-fba is tool to delete 'island cycles' and find unproducible intermediate mets, not for quant growth prediction.
linearStoichInds = sub2ind(size(model.S), nzMat(:,1), nzMat(:,2));
subtractMat = repmat(eps./numNonzeroRowSums, 1, size(model.S,2));

model2=model;
model2.S(linearStoichInds) = model.S(linearStoichInds)-subtractMat(linearStoichInds);