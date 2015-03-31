function [model] = InitModel(model)
    SetParameters;
    model.int_vars = zeros(length(model.rxns),1);
    model.rowlb = zeros(length(model.mets),1);
    model.rowub = zeros(length(model.mets),1);
    model.model_name = 'ecoli';
    model.lb(model.lb==-1000) = MIN_FLUX;
    model.ub(model.ub==1000) = MAX_FLUX;
    
    % exchange reactions added for metabolites for which the MD-FBA constraints
    % cannot be satisfied (i.e. their growth associated demand for de-novo 
    % synthesis)
    %model = addNewExReactions(model);
    
    % Demand reactions represent the new variables d_j in the formulation. One 
    % reaction is added for each metabolite in the model
    model = addDemandReactions(model); 
end