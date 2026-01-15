function GOSPA = MultiSourceSTE_GOSPA(SourceEst,Truth)

    TruthMatrix = [Truth.x' Truth.y' Truth.Q'];
    
    trueNoSource = size(TruthMatrix,1);
    EstNoSource = size(SourceEst,1);
    m = min(trueNoSource, EstNoSource);
    n = max(trueNoSource, EstNoSource);
    
    C = nchoosek(1:n, m);      % all combinations of length m
    permcombination = [];                     % to store all permutations
    
    for i = 1:size(C,1)
        p = perms(C(i,:));      % all permutations of each combination
        permcombination = [permcombination; p];             % append
    end
    
    permNumber = size(permcombination,1);

    if trueNoSource<=EstNoSource
        X = TruthMatrix;
        Y = SourceEst;
    else
        X = SourceEst;
        Y = TruthMatrix;
    end
    
    EstcostVector = zeros(permNumber,1);
    EstPosecostVector = zeros(permNumber,1);

    CuttoffDiswithQ = 50;
    CuttoffDiswithoutQ = 20;
    
    for permIndx = 1:permNumber
        for ns=1:m
            EstcostVector(permIndx) = EstcostVector(permIndx) + min(vecnorm(X(ns,:) - Y(permcombination(permIndx,ns),:),2,2),CuttoffDiswithQ)^2;
            EstPosecostVector(permIndx) = EstPosecostVector(permIndx) + min(vecnorm(X(ns,1:2) - Y(permcombination(permIndx,ns),1:2),2,2),CuttoffDiswithoutQ)^2;
        end
    end
    
    GOSPA.EstError = sqrt(min(EstcostVector) + 0.5*(n-m)*CuttoffDiswithQ^2 );
    GOSPA.EstPoseError = sqrt(min(EstPosecostVector) + 0.5*(n-m)*CuttoffDiswithoutQ^2);
end