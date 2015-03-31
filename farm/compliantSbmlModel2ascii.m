% jmd
% 6.20.12
% fixSBML.m

% maybe can use compList in readCbModel to avoid '__93[_]' ending

% only rxn & met IDs need to be in a restricted character set
function model=compliantSbmlModel2ascii(model)
model.mets=compliantSbmlNames2ascii(model.mets);
model.rxns=compliantSbmlNames2ascii(model.rxns);



