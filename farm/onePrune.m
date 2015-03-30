% jmd
% 9.15.12

% model=limedNC; allowAllTrans=true; extracellSymbol='CCO-EXTRACELLULAR';
% fbaCsense= FBA constrain sense = 'E';
% ran as: v=onePrune(limedNC, true, 'CCO-EXTRACELLULAR');
function v = onePrune(model, allowAllTrans, extracellSymbol, fbaCsense)
if nargin<2
    allowAllTrans=true;
end
if nargin<3
    extracellSymbol='e';
end
if nargin<4
    fbaCsense='E';
end

% if allow all trans, then assume any extracellular metabolite could be
% supplied and removed. this can be implemented by removing the 
% extracellular metabolites from the oneprune optimization. 
% note: this implementation automatically allows flux of purely  
% extracellular rxns, eg invertase, which i think is ok.
if allowAllTrans
    [baseMetNames,compSymbols] = parseMetNames(model.mets);
    extracellMets = model.mets(strcmp(compSymbols, extracellSymbol));
    % 3rd param (false) says don't remove rxns -- want num of rxns to stay
    % constant
    model = removeMetabolites(model,extracellMets,false);
end

nrxns=size(model.S,2);

% lp is struct for solveCobraLP
% opt is: max 0*v + sum t
% subject to: S*v=0; v >= t; 0 <= t <= 1
lp.A=[model.S zeros(size(model.S));
    eye(nrxns) -eye(nrxns)];
lp.b=zeros(1,sum(size(model.S)));
lp.c=[zeros(1,nrxns) ones(1,nrxns)];
lp.lb=zeros(1,2*nrxns);
lp.ub=[inf*ones(1,nrxns) ones(1,nrxns)];
lp.osense=-1;
lp.csense=[repmat(fbaCsense, 1, size(model.S,1)) repmat('G', 1, nrxns)];

sol = solveCobraLP(lp);
v=sol.full(nrxns + (1:nrxns));
