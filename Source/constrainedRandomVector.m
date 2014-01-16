function X = constrainedRandomVector(popsize,LB,UB,A,b)

% LB <= x <= UB
% A*x <= b

LB = LB(:); % make sure column vectors
UB = UB(:);

X = zeros(numel(LB), popsize);

nAttempts = 0;
nFound    = 0;
while nFound < popsize
   
    nAttempts = nAttempts + 1;
    
    x = LB + rand*(UB - LB);
    if all(A*x <= b)
        % feasible
        nFound = nFound + 1;
        X(:,nFound) = x;
    end
    
end


