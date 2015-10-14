% Group Search Optimizer (GSO)
% This function handle the core methodes of GSO
function [bestMember] = Imp_GSO_Func(f)
%     clear
%     clc
    
    functionIndex = f;
    % Unimodal Functions
    if functionIndex == 1, f = @f1; upperBound = -100; lowerBound = 100; maxFitnessEvaluation = 150000; end % Sphere function
    if functionIndex == 2, f = @f2; upperBound = -10; lowerBound = 10; maxFitnessEvaluation = 150000; end   % Schwefel’s Problem 2.22
    if functionIndex == 3, f = @f3; upperBound = -100; lowerBound = 100; maxFitnessEvaluation = 250000; end   % Schwefel’s Problem 1.2
    if functionIndex == 4, f = @f4; upperBound = -100; lowerBound = 100; maxFitnessEvaluation = 150000; end   % Schwefel’s Problem 2.21
    if functionIndex == 5, f = @f5; upperBound = -30; lowerBound = 30; maxFitnessEvaluation = 150000; end   % Generalized Rosenbrock's Function
    if functionIndex == 6, f = @f6; upperBound = -100; lowerBound = 100; maxFitnessEvaluation = 150000; end   % Step Function
    if functionIndex == 7, f = @f7; upperBound = -1.28; lowerBound = 1.28; maxFitnessEvaluation = 150000; end   % Quartic Function i.e. Noise

    % (Complex) Multimodal Functions
    % Multimodal Functions With Many Local Minima
    if functionIndex == 8, f = @f8; upperBound = -500; lowerBound = 500; maxFitnessEvaluation = 150000; end   % Quartic Function i.e. Noise
    if functionIndex == 9, f = @f9; upperBound = -5.12; lowerBound = 5.12; maxFitnessEvaluation = 250000; end   % Generalized Rastrigin’s Function
    if functionIndex == 10, f = @f10; upperBound = -32; lowerBound = 32; maxFitnessEvaluation = 150000; end   % Ackley’s Function
    if functionIndex == 11, f = @f11; upperBound = -600; lowerBound = 600; maxFitnessEvaluation = 150000; end   % Generalized Griewank Function
    if functionIndex == 12, f = @f12; upperBound = -50; lowerBound = 50; maxFitnessEvaluation = 150000; end   % Generalized Penalized Functions 1
    if functionIndex == 13, f = @f13; upperBound = -50; lowerBound = 50; maxFitnessEvaluation = 150000; end   % Generalized Penalized Functions 2
    
    % Group parameters
    groupSize = 48;
    dimensionNumber = 300;
    
    
    % Define parameters
    % f =  fitness function;
    % Upper Bound = upperBound
    % Lower Bound = lowerBound
    % Group Size = groupSize
    % Number of dimensions = dimensionNumber
    
    % Define evolution parameters
    maxIteration      =  10^6; % Maximum number of iterations
    fitnessEvaluation = 0;     % Current fitness evaluation 
    i                 = 1;     % Current iteration number


    % Initilize group position
    upperBound      = repmat(upperBound, groupSize, dimensionNumber);
    lowerBound      = repmat(lowerBound, groupSize, dimensionNumber);
    groupPosition   = lowerBound + ((upperBound - lowerBound) .* rand(groupSize, dimensionNumber));
    groupPosition   = [groupPosition f(groupPosition)];        % Calculate fitness
    upperBound      = repmat(100, 1, dimensionNumber);
    lowerBound      = repmat(-100,1, dimensionNumber);
    bestMember      = min(groupPosition(:,end));    % Best Member
%     fprintf('fitness Evaluation = %u, Best Member = %e\n',fitnessEvaluation,bestMember);

    
    % Initilize group head angle [pi/4 ... pi/4]
    headAngle =  repmat(pi/4, groupSize, dimensionNumber-1);


    % Define GSO key parameters
    producerNumber  = 1;           % Number of Producers
    scroungerNumber = ceil(0.8 * (groupSize - producerNumber));       % Number of Scroungers
    rangerNumber    = groupSize - (producerNumber + scroungerNumber); % Number of Rangers
    a               = round(sqrt(dimensionNumber + 1));               % Turn head back to zero degree
    thetaMax        = pi/a^2;      % Maximum pursuit angle
    aMax            = thetaMax/2;  % Maximum turning angle
    lMax            = sqrt(sum(lowerBound - upperBound)^2);           % Maximum pursuit distance
    
    aIteration      = 0;                                              % Find better area after a iterations
    bestHeadAngle   = zeros(1, dimensionNumber - 1);                  % best head angle

    % Adaptive MAximum maximum pursuit distance
    successCounter = 0; % Success count
    Ps = 0; % Success probability
    wMax = 0.9; % Maximum inertia weight
    wMin = 0.4; % Minimum Inertia weight
    tempPosition = 0;

    while fitnessEvaluation < maxFitnessEvaluation && ...
          i < maxIteration && bestMember(end) ~= 0
        % Select Producers
        if     producerNumber == 1
            [value index]      = min(groupPosition(:,end));
            producerIndex      = index;
%             producer           = groupPosition(producerIndex,:); % Producer members

            % Generate random sequence
            randomSequence     = randperm(groupSize);
            % Omit producer index
            randomSequence     = randomSequence(randomSequence ~= producerIndex);

        elseif producerNumber >= 2   
            [value index]      = sort(groupPosition(:,end));
            producerIndex      = index(1:producerNumber);
%             producer           = groupPosition(producerIndex,:);  % Producer members

            % Generate random sequence
            randomSequence     = randperm(groupSize);
            % Omit producer index
            randomSequence     = randomSequence(randomSequence(~ismember(randomSequence, producerIndex)));
        end

        % Select Scroungers
        scroungerIndex = randomSequence(1:scroungerNumber);
%         scrounger      = groupPosition(scroungerIndex,:);      % Scrounger members

        % Select Rengers
        rangerIndex    = randomSequence(scroungerNumber+1:end);
%         ranger         = groupPosition(rangerIndex,:);         % Ranger members


        % Generate random numbers for producing
        r1 = random('Normal',  0, 1, [producerNumber 1]);                   % Normally  distributed random number
        r2 = random('Uniform', 0, 1, [producerNumber dimensionNumber - 1]); % Uniformly distributed random sequence

        % Perform producing
        for j = 1: producerNumber


            % Randomly sampling three points in the scanning field using
            % equations (2) to (4)

            % Scan at zero degree using equation (2)
            zeroPosition           = zeros(1,dimensionNumber+1);
            zeroPosition(1:end-1)  = groupPosition(producerIndex(j),1:end-1) + ....
                                     r1(j) * lMax * Polar2Cartesian(headAngle(producerIndex(j),:));


            % Restrict search
            % Handle the bounded search space
            outrange                = zeroPosition(1:end-1) > upperBound;
            zeroPosition(outrange)  = groupPosition(producerIndex(j),outrange);
            outrange                = zeroPosition(1:end-1) < lowerBound;
            zeroPosition(outrange)  = groupPosition(producerIndex(j),outrange);
            % Clculate fitness
            zeroPosition(end)       = f(zeroPosition(1:end-1));                          


            % Scan right hand side hypercube using equation (3)
            rightPosition           = zeros(1,dimensionNumber+1);
            rightPosition(1:end-1)  = groupPosition(producerIndex(j),1:end-1) + ...
                                      r1(j) * lMax * Polar2Cartesian(headAngle(producerIndex(j),:) + (r2(j,:) * thetaMax));

            % Restrict search
            % Handle the bounded search space
            outrange                = rightPosition(1:end-1) > upperBound;
            rightPosition(outrange) = groupPosition(producerIndex(j),outrange);
            outrange                = rightPosition(1:end-1) < lowerBound;
            rightPosition(outrange) = groupPosition(producerIndex(j),outrange);
            % Clculate fitness
            rightPosition(end)      = f(rightPosition(1:end-1));


            % Scan left hand side hypercube using equation (4)
            leftPosition            = zeros(1,dimensionNumber+1);
            leftPosition(1:end-1)   = groupPosition(producerIndex(j),1:end-1) + ...
                                      r1(j) * lMax * Polar2Cartesian(headAngle(producerIndex(j),:) - (r2(j,:) * thetaMax));

            % Restrict search
            % Handle the bounded search space
            outrange                = leftPosition(1:end-1) > upperBound;
            leftPosition(outrange)  = groupPosition(producerIndex(j),outrange);
            outrange                = leftPosition(1:end-1) < lowerBound;
            leftPosition(outrange)  = groupPosition(producerIndex(j),outrange);
            % Clculate fitness
            leftPosition(end)       = f(leftPosition(1:end-1));


            % Find the best point with the best resource
            [value index] = min ([zeroPosition(end) rightPosition(end) leftPosition(end) groupPosition(producerIndex(j),end)]);

            % Fly to the best point 
            bestPoint = index;
            switch (bestPoint)
                case {1}
                    % Fly at zero degree
                    groupPosition(producerIndex(j),:) = zeroPosition;
                    headAngle    (producerIndex(j),:) = headAngle(producerIndex(j),:);
                    % Store the last good head angle
                    bestHeadAngle                     = headAngle(producerIndex(j),:);
                    % Restart a counter
                    aIteration                        = 0;
                    % Increment Success Counter
                    successCounter = successCounter + 1;

                case {2}
                    % Fly at right hand side hypercube
                    groupPosition(producerIndex(j),:) = rightPosition;
                    headAngle    (producerIndex(j),:) = headAngle(producerIndex(j),:) + (r2(j,:) * thetaMax);
                    % Store the last good head angle
                    bestHeadAngle                     = headAngle(producerIndex(j),:);
                    % Restart a counter
                    aIteration                        = 0;
                    % Increment Success Counter
                    successCounter = successCounter + 1;


                case {3}
                    % Fly at left hand side hypercube
                    groupPosition(producerIndex(j),:) = leftPosition;
                    headAngle    (producerIndex(j),:) = headAngle(producerIndex(j),:) - (r2(j,:) * thetaMax);
                    % Store the last good head angle
                    bestHeadAngle                     = headAngle(producerIndex(j),:);
                    % Restart a counter
                    aIteration                        = 0;
                    % Increment Success Counter
                    successCounter = successCounter + 1;

                case {4}
                    % Stay at the current Position and turn head angle to a new
                    % randomly generated angel using equation (5)
                    headAngle(producerIndex(j),:)     = headAngle(producerIndex(j),:) + (r2(j,:) *     aMax);
                    % increase a counter by one
                    aIteration                        = aIteration + 1;
                    % if a producer cannot a better area after a iterations, it
                    % will turn its head back to zero degree using equation (6)
                    if aIteration > aMax
                        headAngle(producerIndex(j),:) = bestHeadAngle;
                        % Restart a counter
                        aIteration                    = 0;                    
                    end
            end
        end


            % Generate random numbers for scrounging
            r2 = random('Uniform', 0, 1, [scroungerNumber dimensionNumber - 1]); % Uniformly distributed random sequence
            r3 = random('Uniform', 0, 1, [scroungerNumber dimensionNumber    ]); % Uniformly distributed random sequence

            % Perform scrounging
            for j = 1:scroungerNumber
                
                tempPosition = groupPosition(scroungerIndex(j),end); 
                % Scroungers search for opportunities to join the resources
                % found by the producers
                selectedProducer                         = producerIndex(randi(producerNumber));

                % Scrounger is modeled as a random walk toward the producer
                % using equation (7)
                groupPosition(scroungerIndex(j),1:end-1) = groupPosition(scroungerIndex(j),1:end-1) + ...
                             (r3(j,:) .* (groupPosition(selectedProducer,1:end-1)) - groupPosition(scroungerIndex(j),1:end-1));

                % Restrict search
                % Handle the bounded search space
                outrange                = groupPosition(scroungerIndex(j),1:end-1) > upperBound;
                groupPosition(outrange) = groupPosition(scroungerIndex(j),outrange);
                outrange                = groupPosition(scroungerIndex(j),1:end-1) < lowerBound;
                groupPosition(outrange) = groupPosition(scroungerIndex(j),outrange);
                % Calculate fitness
                groupPosition(scroungerIndex(j),end) = f(groupPosition(scroungerIndex(j),1:end-1)); 


                % The scrounger will keep searching for other opportunities by 
                % turning its head to a new randomly generated angle using (5)
                headAngle(scroungerIndex(j),:)       = headAngle(scroungerIndex(j),:) + (r2(j,:) * aMax);
                
                % Increment Success Counter    
                if tempPosition < groupPosition(scroungerIndex(j),end);
                    successCounter = successCounter + 1;
                end

            end


            % Generate random numbers
            r1 = random('Normal',  0, 1, [rangerNumber 1]);                   % Normally  distributed random number
            r2 = random('Uniform', 0, 1, [rangerNumber dimensionNumber - 1]); % Uniformly distributed random sequence

            % Perform ranging
            for j = 1:rangerNumber       

                tempPosition = groupPosition(scroungerIndex(j),end); 

                % Rangers are disperse and perform random walks to find
                % randomly distributed resources
                % Rangers generate a random head angle using equation (5)
                headAngle(rangerIndex(j),:) = headAngle(rangerIndex(j),:) + (r2(j,:) * aMax);

                % Rangers choosse a random distance using equation (8)
                randomDistance = a * r1(j) * lMax;

                % Rangers move to the new point using equation (9)
                groupPosition(rangerIndex(j),1:end-1) = groupPosition(rangerIndex(j),1:end-1) + ...
                                randomDistance * Polar2Cartesian(headAngle(rangerIndex(j),:));


                % Restrict search
                % Handle the bounded search space
                outrange                = groupPosition(rangerIndex(j),1:end-1) > upperBound;
                groupPosition(outrange) = groupPosition(rangerIndex(j),outrange);
                outrange                = groupPosition(rangerIndex(j),1:end-1) < lowerBound;
                groupPosition(outrange) = groupPosition(rangerIndex(j),outrange);
                % Calculate fitness
                groupPosition(rangerIndex(j),end) = f(groupPosition(rangerIndex(j),1:end-1));
                
                % Increment Success Counter    
                if tempPosition < groupPosition(scroungerIndex(j),end);
                    successCounter = successCounter + 1;
                end
            end
            
            
            % Calculate success percentage
            Ps = successCounter / groupSize;
            lMax = (wMax - wMin) * Ps  + wMin;
            successCounter = 0;
            

        % Calculate best producer and show its fitness
        bestMember = [bestMember min(groupPosition(:,end))];
 
        % Calculate number of fitness evaluation
        fitnessEvaluation = fitnessEvaluation + groupSize;
        % insrease iteration number
        i = i + 1;
        
        if mod(i,100) == 1
            fprintf('fistness Evaluation = %u, Best Member = %e, lMax = %e\n',fitnessEvaluation,bestMember(i), lMax);
        end
        
        
        if i > 100
            if bestMember(i - 50) == bestMember(i)
                break;
            end
        end
%         
%         if bestMember(end) == 0
%             break;
%         end
    end
end