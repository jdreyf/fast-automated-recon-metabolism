function [new_model]= addNewExReactions(model)
    
    new_model = model;
    % Add input exchange reactions for the following proteins / trna!:
    % alpp                           (applipoprotein)
    % apoACP                      (apoprotein [acyl carrier protein])
    % trdox                        (Oxidized thioredoxin)
    % trnaglu                   (tRNA (Glu))
    % fldrd[c] 761 Flavodoxin (reduced)[c]
    % grxrd[c] 922 Glutaredoxin (reduced)[c]
    % 551 D-Carnitine[c]

    protein1 = strmatch('alpp',model.mets);
    protein2 = strmatch('apoACP',model.mets);
    protein3 = strmatch('trdox',model.mets);
    trna = strmatch('trnaglu',model.mets);
    mets = [protein1;protein2;protein3;trna;761;922;551];
    cols = [1:length(mets)];
    vals = -1*ones(length(mets),1);
    mat = sparse(mets,cols,vals,length(model.mets),length(mets));
    new_model.S = [model.S mat];
    new_model.lb = [model.lb; -8*ones(length(mets),1)];
    new_model.ub = [model.ub; zeros(length(mets),1)];
    new_model.c = [model.c; zeros(length(mets),1)];
    new_model.int_vars = [model.int_vars; zeros(length(mets),1)];
    New_rxnGeneMat = sparse (length(mets),length(model.genes));
    new_model.rxnGeneMat = [model.rxnGeneMat; New_rxnGeneMat];
    index = length(model.rxns);
    new_model.rxns{index+1} = 'New_ex_alpp';
    new_model.rxns{index+2} = 'New_ex_apoACP';
    new_model.rxns{index+3} = 'New_ex_trdox';
    new_model.rxns{index+4} = 'New_ex_trnaglu';
    new_model.rxns{index+5} = 'New_ex_fldrd';
    new_model.rxns{index+6} = 'New_ex_grxrd';
    new_model.rxns{index+7} = 'New_ex_crn-DASH-D[c]';
end