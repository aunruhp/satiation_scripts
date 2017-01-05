
function [bslshuff] = randshuff_hcoeff(nstims,rspdur,rawdata,mfr,fast,nshuff)
%
% randshuff is a function that takes all spikes from within one whole
% recording session (rawdata) and from it creates nshuff (i.e. usually 1000
% due to time constraints) randomly generated responses irrespective of whether
% there actually was a real response during some part or not. To each one
% of these shuffled responses randshuff applies adptive kernel
% smoothing (ssvkernel) and then calculates a master response M across all
% of them, from which one can then later calculate the h-coefficient.
% To execute ssvkernel 1000 times takes a very long time on a normal PC, so
% randshuff also offers a 'fast' method (fast = 1). Instead of cutting out
% 1000 segments of data of a certain length and then applying ssvkernel to
% each segment, the fast method cuts out 100 segments of 10 times that
% length. It then only needs to apply ssvkernel 100 times instead of 1000
% times. After running ssvkernel 10 times randshuff then cuts each smoothed
% segment into 10 parts, which again results in 1000 segments. Running
% randshuff in fast mode is not as good as running it in 'slow' (fast = 0)
% mode but in most cases it's ok as long as your interstimulus intervalls
% are at least 10 times longer then your baseline periods.
%
% Input:
%
% nstims:
% nstims is the number of trials or stimuli that were run during the whole
% experiments
%
% rspdur:
% rspdur is the response duration in ms, i.e. the duration of the individual
% segments that will make up your bootstraped responses
%
% rawdata
% rawdata is a vector containng timestamps for each spike that was recorded
% during yoru experiment.
%
% mfr
% mfr is the mean firing rate for this neuron across the whole
% experiment
%
% fast
% fast = 1 runs ranshuff in fast mode. fast = 0 runs randshuff in slow mode

dbstop if error
if nargin < 5
    fast = 1;
end
if nargin < 6
    nshuff = 1000;
end
kernshuff = randperm(nshuff);
shuffsy = zeros(nstims,1) * NaN;
nsegs = floor(((floor(rawdata(end))-2000)-(ceil(rawdata(1,1))))/(rspdur*10));   % divide whole time of raw data in 10 segments, nseg is number of segments   % Add 2 x 1000 = 2000 ms here so that you can always subtract them again afterwards to get rid of any edge artifact in the ssvkernel
if fast
    for i = 1:nshuff/10
        clc
        disp([num2str(((i-1)*10)+1),' to ',num2str(i*10),' out of ',num2str(nshuff),' shuffles being processed with the "fast" method...'])
        permind1 = randperm(nsegs*10);
        permind2 = randperm(nsegs); %random numbers for starting point of each segment
        for j = 1:nstims
            shuffsy(j,1) = length(rawdata((rawdata > ((permind1(1,j)-1)*rspdur) & (rawdata <= (((permind1(1,j)-1)*rspdur)+rspdur)))));
            if j == 1 % Prep for the adaptive kernel analysis
                forkern = (rawdata((rawdata > ((permind2(1,j)-1)*(rspdur*10))) & (rawdata <= ((permind2(1,j)-1)*(rspdur*10))+(rspdur*10)+2000)) - ((permind2(1,j)-1)*(rspdur*10)));
            else  % Only do this 100 times, but with a x10 longer duration to save time (instead of 100 times)
                forkern = cat(2,forkern,(rawdata((rawdata > ((permind2(1,j)-1)*(rspdur*10))) & (rawdata <= ((permind2(1,j)-1)*(rspdur*10))+(rspdur*10)+2000)) - ((permind2(1,j)-1)*(rspdur*10))));
            end
        end
        [prekern{i}(2,:),prekern{i}(1,:)] = ssvkernel(forkern.*10^(-3),(0.01:0.01:(((rspdur*10)+2000)/1000))); % Execute your adaptive kernel analysis
        prekern{i}(2,:) = prekern{i}(2,:) .* 10 .* (squeeze(nanmean(shuffsy))*(1000/rspdur));
        for j = 1:10
            k = ((i-1)*10)+j;
            bslshuff.kern{kernshuff(k)}(1,:) = (0.01:0.01:rspdur/1000);
            bslshuff.kern{kernshuff(k)}(2,:) = prekern{i}(2,(j-1)*(rspdur/10)+1:j*(rspdur/10));      % Reverse the normaliztion applied by ssvkernel and transform to Hz. Also cut off the first and last 50 ms, as these have edge artifact. You can do so because above you've added 100 ms to the baseline, so youre going to end up with the right length again.
        end
    end
else
    for i = 1:nshuff
        clc
        disp([num2str(i),' out of ',num2str(nshuff),' shuffles being processed with the "slow" method (this might take a while)...'])
        permind = randperm(nsegs*10); %10 more segments, 10 times shorter
        for j = 1:nstims
            shuffsy(j,1) = length(rawdata((rawdata > ((permind(1,j)-1)*rspdur) & (rawdata <= (((permind(1,j)-1)*rspdur)+rspdur)))));
            if j == 1   % Prep for the adaptive kernel analysis
                forkern = (rawdata((rawdata > ((permind(1,j)-1)*rspdur)) & (rawdata <= ((permind(1,j)-1)*rspdur)+rspdur+2000)) - ((permind(1,j)-1)*rspdur))';
            else
                forkern = cat(2,forkern,(rawdata((rawdata > ((permind(1,j)-1)*rspdur)) & (rawdata <= ((permind(1,j)-1)*rspdur)+rspdur+2000)) - ((permind(1,j)-1)*rspdur))');
            end
        end
        [prekern{i}(2,:),prekern{i}(1,:)] = ssvkernel(forkern.*10^(-3),(0.01:0.01:((rspdur+2000)/1000))); % Execute your adaptive kernel analysis
        prekern{i}(2,:) = prekern{i}(2,:) .* 10 .* (squeeze(nanmean(shuffsy))*(1000/rspdur)); %in Hz
        bslshuff.kern{kernshuff(i)}(1,:) = (0.01:0.01:rspdur/1000);
        bslshuff.kern{kernshuff(i)}(2,:) = prekern{i}(2,:);
    end
end
clc

%% Now go through your shuffled adaptively smoothed baselines/responses and look for responses
bslshuffmaxes = zeros(10001,nshuff);
bslshuffminses = bslshuffmaxes;
for a = 1:nshuff  % for each of your shuffled samples
    kern = (bslshuff.kern{a}(2,:)./mfr)-1;                                  % Normalize it to your channel's mean firing rate and subtract 1 so that your baseline value is = 0
    indmax = find(kern == max(kern),1,'first');                             % Find the areas for your local maximum and save them into the maxes variable.
    if ~isempty(indmax) && kern(1,indmax) > 0                               % just to make sure your maximum is above the mean firing rate of this channel, as else this analysis makes no sense
        firstdec = norm(floor(norm(kern(1,indmax))*10)/10);                 % Find out what the next decimal point is between your max and the channel baseline (in units of channel baseline firing rate, e.g. if the max firing rate is 1.568 times the mean firing rate of the channel, then your next decimal point towards that mean firing rate is 1.5)
        u = indmax(1); v = u;
        for b = 1:round(firstdec*10)+1
            while u > 1 && kern(1,u) > roundn(firstdec-((b-1)*0.1),-1)      % This travels along the curve towards the left (heading downwards from the local maximum) until it hits the next decimal multiple of the mean firing rate of the whole channel
                u = u-1;
            end
            while v < length(kern) && kern(1,v) > roundn(firstdec-((b-1)*0.1),-1)   % This travels along the curve towards the right (heading downwards from the local maximum) until it hits the next decimal multiple of the mean firing rate of the whole channel
                v = v+1;
            end
            bslshuffmaxes(end-(round(firstdec*10)+1)+b,a) = trapz((kern(1,u:v))-roundn(firstdec-((b-1)*0.1),-1)).*0.01; % "-firstdec" makes sure the only the area between the next decimal multiple of the baseline firing rate and the peack that reaches just above this threshold is calculated. The factor 0.1 is there because the spacing of our x-axis is 0.1, see help trapz for more on this
        end
    end
    indmin = find(kern == min(kern),1,'first');                             % Find the areas for your local minimum and save them into the minses variable.
    if ~isempty(indmax) && kern(1,indmin) < 0                               % just to make sure your maximum is above the mean firing rate of this channel, as else this analysis makes no sense
        firstdec = norm(floor(norm(kern(1,indmin))*10)/10);                 % Find out what the next decimal point is between your max and the channel baseline (in units of channel baseline firing rate, e.g. if the max firing rate is 1.568 times the mean firing rate of the channel, then your next decimal point towards that mean firing rate is 1.5)
        u = indmin(1); v = u;
        for b = 1:round(firstdec*10)+1
            while u > 1 && kern(1,u) < -roundn(firstdec-((b-1)*0.1),-1)     % This travels along the curve towards the left (heading upwards from the local minimum) until it hits the next decimal multiple of the mean firing rate of the whole channel
                u = u-1;
            end
            while v < length(kern) && kern(1,v) < -roundn(firstdec-((b-1)*0.1),-1) % This travels along the curve towards the right (heading upwards from the local minimum) until it hits the next decimal multiple of the mean firing rate of the whole channel
                v = v+1;
            end
            bslshuffminses(end-(round(firstdec*10)+1)+b,a) = norm(trapz(kern(1,u:v) + roundn(firstdec-((b-1)*0.1),-1)).*0.01); % "+(firstdec-((b-1)*0.1))" makes sure that only the area between the next decimal multiple of the baseline firing rate and the peack that reaches just below this threshold is calculated. The factor 0.1 is there because the spacing of our x-axis is 0.1, see help trapz for more on this
        end
    end
end
bslshuffmaxes = bslshuffmaxes(2:end,:) - bslshuffmaxes(1:end-1,:);             % This subtracts the row above from each row. Before doing so your area measurements represent the area above each decima multiple of the baseline firing rate. After doing this, your area measurements are now corrected to reflect only the area inbetween this and the next following decimal. In other words, instead of reporting how large the area is above e.g. 1.7 * meanfr(channel), we now are reporting how large the area is inbetween 1.7 * meanfr(channel) and 1.8 * meanfr(channel), i.e. the maximum in each band
bslshuffminses = bslshuffminses(2:end,:) - bslshuffminses(1:end-1,:);
bslshuff.maxes = max(bslshuffmaxes,[],2);
bslshuff.minses = max(bslshuffminses,[],2);
end