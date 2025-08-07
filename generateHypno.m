function sleepStages = generateHypno(totalSleep, remProportion)
    % Define parameters
    cycleDuration = 90*30; % Each sleep cycle is 90 minutes
    numCycles = floor(totalSleep / cycleDuration); % Number of full cycles
    if numCycles < 1
        error('Sleep duration must be at least one full cycle (90 minutes).');
    end
    
    % Define REM proportion increasing across cycles
    baseREM = repelem(remProportion,numCycles); % Rough baseline REM per cycle
    remmodulation = linspace(0.5, 1.5, numCycles); % Increase REM in later cycles
    remFractions = baseREM .* remmodulation;
    % Construct sleep vector
    sleepStages = [];
    
    for i = 1:numCycles
        cycle = zeros(1, cycleDuration); % Start with all NREM
        remDuration = round(remFractions(i) * cycleDuration)/100; % Compute REM duration
        cycle(end-remDuration+1:end) = 2; % Assign REM to end of cycle
        sleepStages = [sleepStages, cycle]; % Append to final vector
    end
    sleepStages(sleepStages==0) = 3;
end