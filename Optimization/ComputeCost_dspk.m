% Currents that modulate spike shape dynamics on multiple timescales induce
% ramping bursts in model respiratory neurons
% Code written by : Dr. Victor Matveev, NJIT

function Cost = ComputeCost_dspk( Tlist, Vlist, plotFlag )

%    Inputs: "Tlist": time array; "Vlist": V(t) array
%   * Process spike times, detect bursts, and 
%   * Return "Cost" that favors ramp-up bursts
%   * Make plots unless "plotFlag" is zero

minV  = min(Vlist);               % Maximal voltage
maxV  = max(Vlist);               % Minimal voltage
rng    = 2 : numel(Vlist)-1;      % Indices of V(t) array, except endpoints

% Note a trick: "TestMinMax" positive only at local min or max of V(t):
TestMinMax = (Vlist(rng) - Vlist(rng-1)) .* (Vlist(rng) - Vlist(rng+1));   
iMaxMin    = find(TestMinMax > 0);        

TestSpike  = Vlist(iMaxMin) - Vlist(iMaxMin-1);
spikeIndex = iMaxMin(find(TestSpike > 0))+1;
spikeTimes = Tlist(spikeIndex); 
spikeLevels= Vlist(spikeIndex);

TestTrough  = Vlist(iMaxMin-1)-Vlist(iMaxMin);
minIndex = iMaxMin(find(TestTrough > 0));
minTimes = Tlist(minIndex); 
minLevels = Vlist(minIndex);

if plotFlag               % Put circles at spike and trough times: 
   hold off; plot(Tlist, Vlist); hold on;
   plot(spikeTimes, spikeLevels,'mo');
   plot(minTimes, minLevels,'bo');
   drawnow;
end

ISIs     = diff( spikeTimes );    % Compute inter-spike intervals
maxISI   = max(ISIs);             % Compute maximal ISI
minISI   = min(ISIs);             % Compute minimal ISI  
spikeISI = 0.5*(minISI + maxISI); % Cut-off btwn burst and inter-burst ISI

if maxISI / minISI < 3            % Require that maxISI > 3 * minISI
    Cost = 100;                   % Otherwise set cost to max and return
    title('No bursts detected'); 
    fprintf('No burst detected\n');
    drawnow;
    return;
end

burstStart = [];                  % Initialize burst starting ISI array 
burstEnd   = [];                  % Initialize burst ending ISI array 
inBurst    = false;               % True if ISI belongs to a burst

for k = 1 : numel( ISIs )         % Test the lengths of all ISIs
    if ISIs(k) < spikeISI         % ISI small enough to belong to a burst
        if ~inBurst               % If we are not in a burst yet, then:
            burstStart = [burstStart, k];  % Store first ISI index
            inBurst    = true;             % Set in-burst flag to "true"
        end
    else                          % Else ISI is a larger inter-burst ISI
        if inBurst                % If we were in a burst until now, then:
            burstEnd = [burstEnd, k];   % Store last ISI index
            inBurst  = false;           % Reset in-burst flag
        else                      % Otherwise, this is a second long ISI, so:
            Cost = 100;           % Set cost to maximal value, and return
            title(' Not a clean burst '); 
            fprintf(' Not a clean burst \n');
            drawnow; return;
        end
    end
end

if numel(burstEnd) < 3                      % Need there to be at least three bursts 
    Cost = 0;                               % Zero cost will trigger extra computation
    title(' Not enough bursts: extend time interval');
    fprintf(' Not enough bursts: extend time interval\n');
    drawnow; return;
end

if(numel(burstStart(2):burstEnd(2))~=numel(burstStart(3):burstEnd(3)))
    Cost = 100;           % Set cost to maximal value, and return
    title(' Not a clean burst '); 
    fprintf(' Not a clean burst \n');
    drawnow; return;
end

spikeFirst = burstStart(2);                 % Index of first spike in the 2nd burst
spikeLast  = burstEnd(2);                   % Index of last  spike in the 2nd burst
rng        = spikeFirst : (spikeLast-1);    % Range of ISIs within 2nd burst

if(numel(rng)<5)
    Cost=100;
    title(' Not enough spikes in burst');
    fprintf(' Not enough spikes in burst \n');
    drawnow; return;
end

if(minLevels(rng)<(spikeLevels(rng)-1))
else 
    Cost=100;
    title('Depolarization block in burst');
    fprintf('Depolarization block in burst \n');
    drawnow; return;
end

if plotFlag                                 % Put black lines around burst #2
    plot(spikeTimes( spikeFirst ) * [1 1], [minV maxV], 'k-');
    plot(spikeTimes( spikeLast  ) * [1 1], [minV maxV], 'k-');
    drawnow;
end

Ratio1 = ISIs(rng) ./ ISIs(rng-1);                          % Ratio of successive ISIs 
Ratio2 = minLevels(rng)./ minLevels(rng-1);                 % Ratio of successive plateaus 
cost1 = mean(Ratio1) + 4*sum(Ratio1(end-6:end));            % Mean ratio of successive ISIs
cost2 = mean(Ratio2) + 4*sum(Ratio2(end-6:end)) ;           % Mean ratio if successive V-troughs
cost3 = 0.1 * max(ISIs(rng)) / maxISI;                      % Penalize (max-ISI/max-burst-ISI) ratio
cost4 = 2*min(spikeLevels(rng))./ spikeLevels(spikeFirst);  % Penalize spike peak reduction
Cost  = sqrt(cost1^2 + cost2^2 + cost3^2 + cost4^2);  % Combine above into a single Cost

str = sprintf(' Cost: [%g %g %g %g] => %g \n', cost1, cost2, cost3, cost4, Cost);
fprintf('%s', str);

if plotFlag
    title(str);
    drawnow;
end

end  % %%%%%%%%%%%%%%%%%% END OF FUNCTION %%%%%%%%%%%%%%%%%%%%

