% could accept transRxns to subset model
function [pureTransRxns,nonPureTransRxns] = findPureTransRxns(model, inclExc)
%findPureTransRxns identifies all reactions in a model that only transport substrates but do not
% change them in any way. E.g. ATP synthase involves multiple compartments, but does more than
% simply transport, so it's not a 'pure' transport reaction.
%
% [transRxns,nonTransRxns] = findTransRxns(model,inclExc)
%
%INPUT
% model             COBRA model structure
%
%OPTIONAL INPUT
% inclExc           include exchange reactions as transport?
%                   (Default = false)
%
%OUTPUT
% pureTransRxns         all pure transport reactions in model
% nonPureTransRxns      all not pure transport reactions in model
%
% Jonathan Dreyfuss, June 2012

if nargin < 2
    inclExc = false;
end

if inclExc
    % findExcRxns returns boolean vector
    isExc = findExcRxns(model,inclObjFlag,irrevFlag);
else
    isExc=zeros(1, length(model.rxns));
end

isNonexchPureTrans = zeros(1,length(model.rxns));

baseMetNames=arrayfun(@parseMetNames, model.mets);

for i = 1:length(model.rxns)
    nonzeroInd = find(model.S(:,i));
	stoichCoeffs=model.S(nonzeroInd,i);
	baseMetsTmp=baseMetNames(nonzeroInd);
    uniqueBaseMetsTmp=unique(baseMetsTmp);
    
    % reaction is pure transport if and only if the sum of coeffs 
    % corresponding to each of its baseMetNames are 0
	baseMetCoeffs=zeros(1,length(uniqueBaseMetsTmp));
    for j=1:length(uniqueBaseMetsTmp)
		metStoichCoeffs=stoichCoeffs(ismember(baseMetsTmp,uniqueBaseMetsTmp{j}));
        baseMetCoeffs(j)=sum(metStoichCoeffs);
    end
	
    if all(baseMetCoeffs==0)
        isNonexchPureTrans(i)=1;
    else
        isNonexchPureTrans(i)=0;
    end
end

pureTransRxns = model.rxns(isNonexchPureTrans==1 | isExc==1);
nonPureTransRxns = setdiff(model.rxns, pureTransRxns);
