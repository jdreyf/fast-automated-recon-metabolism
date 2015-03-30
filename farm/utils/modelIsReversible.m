% jmd
% 9.12.12
% isModelReversible.m

% not all models have a model.reversibleModel slot, & when it's there,
% it need not be accurate
% assumes ub>=lb for all rxns
function revModel = modelIsReversible(model)
% model is reversible if and only if there exists a rxn whose ub>0 & lb<0
revModel=any(sign(model.lb).*sign(model.ub)<0);