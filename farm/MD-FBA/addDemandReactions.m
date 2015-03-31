function [model] = addDemandReactions (model)
    SetParameters;
    diag_vec = (-1) * ones(length(model.mets),1);
    mat_ex = sparse(diag(diag_vec));
    model.S = [model.S mat_ex];
    model.lb = [model.lb; zeros(length(model.mets),1)];
    model.ub = [model.ub; MAX_FLUX*ones(length(model.mets),1)];
    model.c = [model.c; zeros(length(model.mets),1)];
    model.int_vars = [model.int_vars;zeros(length(model.mets),1)];
    model.rxns = [model.rxns; cellstr(num2str((1:length(model.mets))'))];
end