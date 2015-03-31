function [milp] = CreateMDFBA (model)

    irrelevant_mets = [];
    name = 'CBFBA';
    print_level = 0;

    SetParameters;
    milp = model;  
    mets_len = size(model.S,1);
    demand_rxns_len = mets_len;
    orig_rxns_len = size(model.S,2) - demand_rxns_len;
    orig_S = model.S(:,1:orig_rxns_len);
    mets_indexes = (1:mets_len)';
    relevant_mets = setdiff(mets_indexes,irrelevant_mets);
    int_variables_amount = length(relevant_mets);
    
    % add new variables
    % zero matrix for the new integer vars
    mat_var = sparse(mets_len, int_variables_amount);
    milp.S = [model.S mat_var];
    
    % add new ex reactions constraints
    % zero matrix below original S
    mat1 = sparse(int_variables_amount, orig_rxns_len);
    % constraint coeffs matrix for the demand reactions flux variables
    index_vec = (1:int_variables_amount)';
    mat2 = sparse(index_vec, relevant_mets, ones(int_variables_amount,1),...
                                    int_variables_amount, mets_len);
    % constraint coeffs matrix for the new integer variables
    diag_vec = -EPSILON_ACTIVE * ones(int_variables_amount,1);
    mat3 = sparse(diag (diag_vec));
    milp.S = [milp.S; mat1 mat2 mat3];
    
    % add inner reactions constraints
    for i=1:int_variables_amount
        % find the reactions which involve relevant metabolite i
        rxns_involved = find(orig_S(relevant_mets(i),:)~=0);
        rxns_involved = rxns_involved';
        rxns_inv_len = length(rxns_involved);
        % add positive constraint
        index_vec = (1:rxns_inv_len)';
        mat1 = sparse(index_vec, rxns_involved, ones(rxns_inv_len,1),...
                                        rxns_inv_len, orig_rxns_len);
        mat2 = sparse(rxns_inv_len, mets_len);
        vec = i*ones(rxns_inv_len,1);
        mat3 = sparse(index_vec, vec, INT_ACTIVITY*ones(rxns_inv_len,1),...
                                        rxns_inv_len, int_variables_amount);
        milp.S = [milp.S; mat1 mat2 mat3];
        % add negative constraint
        index_vec = (1:rxns_inv_len)';
        values = (-1)*ones(rxns_inv_len,1);
        mat1 = sparse(index_vec, rxns_involved, values,...
                                        rxns_inv_len, orig_rxns_len);
        mat2 = sparse(rxns_inv_len, mets_len);       
        vec = i*ones(rxns_inv_len,1);
        mat3 = sparse(index_vec, vec, INT_ACTIVITY*ones(rxns_inv_len,1),...
                                        rxns_inv_len, int_variables_amount);
        milp.S = [milp.S; mat1 mat2 mat3];
    end    

    % update lb and ub for new variables
    milp.lb = [milp.lb; zeros(int_variables_amount,1)];   
    milp.ub = [milp.ub; ones(int_variables_amount,1) ];
   
    % update rowlb and rowub for new variables
    % ub and lb for the rows of inner reaction constraints 
    inner_rxn_constraints_len = size(milp.S,1) - mets_len - int_variables_amount;
    % row lb
    vec = (-1)*EPSILON_ACTIVE*ones(inner_rxn_constraints_len,1);
    milp.rowlb = [model.rowlb; ... 
                  zeros(int_variables_amount,1); ... % demand rxns constraints
                  vec]; % Inner rxns constraints
    % row ub
    vec1 = MAX_FLUX*ones(int_variables_amount,1);      % demand rxns constraints
    % have to be 2*MAX_FLUX in order to maintain positive flux in the inner
    % reactions
    vec2 = 2*MAX_FLUX*ones(inner_rxn_constraints_len,1); % inner rxns constraints
    milp.rowub = [model.rowub; vec1;vec2];
    
    % update int_vars
    milp.int_vars = [milp.int_vars; ones(int_variables_amount,1)];
    % update c vector
    milp.c = [milp.c; zeros(int_variables_amount,1)];  
    % update model name
    milp.model_name = name;
    if (print_level > 0)
        save ([name '_milp.mat'], '-struct','milp');
    end
end
