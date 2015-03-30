% jmd
% 9.11.12
% findNondiluteRxns.m

function nonDiluteRxns = findNondiluteRxns(model, allTrans)
% get names of rxns that shouldn't be diluted by limed-FBA
if nargin<2
    allTrans=false;
end

% need irrev model
if modelIsReversible(model)
    error('limedFBA:rev', 'an irreversible model is required');
end

if allTrans
    trans=findTransRxns3(model);
else
    trans=findPureTransRxns(model);
end

isTrans=ismember(model.rxns, trans);

% this includes biomass rxn, which has model.c~=0
isExc = findExcRxns(model,true,true); % CHECK

nonDiluteRxns=model.rxns(isExc | isTrans);