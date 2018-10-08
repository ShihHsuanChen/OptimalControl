% Optimizer.m
% Nelder-Mead Method

%function [para_best, record_x, record_y] = Optimizer(costFunction,parameters)
function [para_best,cost_best] = Optimizer(costfunction,parameters)
% costfunction : a function to calculate infidelity with a parameter set.
%                Should be a object of function such as @Myfun
% parameters   : parametersets [1xn array]
% para_best    : optimized result (parameter)
% cost_best    : optimized result (cost)

    X0 = reshape(parameters,1,[]);
    costfunction = fcnchk(costfunction);

    dim = length(X0); % dimension of the problem

    X_arr = zeros(dim+1,dim);
    F_arr = zeros(dim+1,1);
    % set up adaptive parameters
    alpha = 1;
    beta  = 1+2/dim;
    gamma = 0.75-0.5/dim;
    delta = 1-1/dim;
    max_feval = 100000;
    tol = 10^-6;
    
    % Construct the initial simplex: Large initial simplex is used.
    scalefactor = min(max(max(abs(X0)),1),10);

    D0 = eye(dim);
    D0(dim+1,:) = (1-sqrt(dim+1))/dim*ones(1,dim);
    iterNo = -1;
    for i = 1:dim+1
        X = X0 + scalefactor * D0(i,:);
        X_arr(i,:) = X;
        F_arr(i) = feval(costfunction,X);
    end

    ct = dim+1;

    % Create start point
    [F_arr,Index] = sort(F_arr);
    X_arr = X_arr(Index,:);

    Cost_Container = [];
    iterNo = 0;
    
    % Main iteration
    while max(max(abs(X_arr(2:dim+1,:)-X_arr(1:dim,:)))) >= scalefactor*tol 
        if ct > max_feval
            break;
        end
        
        M = mean(X_arr(1:dim,:));  % Centroid of the dim best vertices
        Xref = (1+alpha)*M - alpha * X_arr(dim+1,:);
        X = Xref;
        Fref = feval(costfunction,X);
        ct = ct + 1;
        
        if Fref < F_arr(1)
            % expansion
            X = (1+alpha*beta)* - alpha*beta*X_arr(dim+1,:);
            Fexp = feval(costfunction,X);
            ct = ct+1;
            if Fexp < Fref
                X_arr(dim+1,:) = X;
                F_arr(dim+1)   = Fexp;
            else
                X_arr(dim+1,:) = Xref;
                F_arr(dim+1)   = Fref;
            end
            
        elseif Fref < F_arr(dim)
            % accept reflection point
            X_arr(dim+1,:) = Xref;
            F_arr(dim+1)   = Fref;
            
        elseif Fref < F_arr(dim+1)
            % Outside contraction
            Xoc   = (1+alpha*gamma)*M-alpha*gamma*X_arr(dim+1,:);
            X_arr = Xoc;
            Foc   = feval(costfunction,X);
            ct=ct+1;           
            if Foc <= Fref
                X_arr(dim+1,:) = Xoc;
                F_arr(dim+1)   = Foc;
            else
                % shrink
                for i=2:dim+1
                    X = X_arr(1,:)+ delta*(X_arr(i,:)-X_arr(1,:));
                    X_arr(i,:) = X;
                    F_arr(i)   = feval(costfunction,X);
                end
                    ct = ct+dim;
            end
        else
            %inside contraction
            Xic=(1-gamma)*M+gamma*X_arr(dim+1,:);
            X = Xic;
            Fic = feval(costfunction,X);
            ct = ct+1;
            if Fic<F_arr(dim+1)
                X_arr(dim+1,:) = Xic;
                F_arr(dim+1)   = Fic;
            else
                % shrink
                for i=2:dim+1
                    X = X_arr(1,:)+ delta*(X_arr(i,:)-X_arr(1,:));
                    X_arr(i,:) = X;
                    F_arr(i)   = feval(costfunction,X);
                end
                ct = ct+dim;
            end
        end
        [F_arr,Index] = sort(F_arr);
        fprintf('%d\t%f\n',iterNo,F_arr(1));
        iterNo = iterNo + 1;
    end
    
    [F_arr,Index] = sort(F_arr);
    
    Cost_Container = [Cost_Container F_arr(1)];
    plot(Cost_Container,'o-')
    
    para_best = X_arr(Index,:);
    para_best = para_best(1,:);
    cost_best = F_arr(1);
end