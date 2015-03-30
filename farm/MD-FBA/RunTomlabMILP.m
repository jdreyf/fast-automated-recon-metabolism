function [milp_model] = RunTomlabMILP(milp_model, print_level, tilim)
    if nargin<2
        print_level=0;
    end
    if nargin<3
        tilim=400;
    end
    % Creates and solves mixed-integer linear programming problems using the
    % TOMLAB format from the standard form of milp problem of a metabolic model
    %Name=milp_model.name;
    Name =milp_model.model_name;
    % Problem formulated as a maximum problem
    A = milp_model.S;
    b_U = milp_model.rowub;
    b_L = milp_model.rowlb;
    c = milp_model.c;
    c=-c;
    [m,n] = size(A);
    x_L = milp_model.lb;
    %%x_L(find(x_L == -1000)) = -1.0E+20;
    x_U = milp_model.ub;
    %%x_U(find(x_U == 1000)) = 1.0E+20;
    x_0 = zeros(n,1);
    if (print_level > 0)
        fprintf('MILP problem. Variables %d. Knapsacks %d\n',n,m);
    end
    IntVars = milp_model.int_vars; 
    x_min = x_L; x_max = x_U; f_Low = -1E7; % f_Low <= f_optimal must hold
    f_opt = [];%Optimal function value(s), if known (Stationary points) not used
    nProblem = []; % Problem number not used
    fIP = []; % Do not use any prior knowledge
    xIP = []; % Do not use any prior knowledge
    setupFile = []; % Just define the Prob structure, not any permanent setup file
    x_opt = []; % The optimal integer solution is not known
    VarWeight = []; % No variable priorities, largest fractional part will be used
    KNAPSACK = 0; % Run without the knapsack heuristic
    % Assign routine for defining a MIP problem.
    Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name, setupFile, ...
    nProblem, IntVars, VarWeight, KNAPSACK, fIP, xIP, ...
    f_Low, x_min, x_max, f_opt, x_opt);
    Prob.optParam.IterPrint = 0; % Set to 1 to see iterations.
    Prob.Solver.Alg = 2; % Depth First, then Breadth search
    % time out
    Prob.MIP.cpxControl.TILIM = tilim;
    % display
    Prob.MIP.cpxControl.TUNINGDISPLAY = 3;

    % Calling driver routine tomRun to run the solver.
    % The 1 sets the print level after optimization.
    milp_model.tomlab_result = tomRun('cplex', Prob, print_level);
    milp_model.result_vector = milp_model.tomlab_result.x_k;
    milp_model.result_opt = -milp_model.tomlab_result.f_k;
    milp_model.result_status = milp_model.tomlab_result.ExitFlag;
    milp_model.result_status_text = milp_model.tomlab_result.ExitText;
    if print_level > 0
        fprintf('Opt val: %d\nExit flag: %d\nExit text: %s\n',milp_model.result_opt,milp_model.result_status,milp_model.result_status_text);
    end
end