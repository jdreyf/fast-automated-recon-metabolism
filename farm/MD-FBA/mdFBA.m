function [milp_res] = mdFBA(model, print_level, tilim)
    if nargin<2
        print_level=0;
    end
    if nargin<3
        tilim=400;
    end
    % initialize model
    model = InitModel(model);
    
    % create the MD-FBA problem
    milp = CreateMDFBA(model);
    
    % solve a MD-FBA MILP problem using the CPLEX solver via Tomlab
    milp_res = RunTomlabMILP(milp, print_level, tilim);
end