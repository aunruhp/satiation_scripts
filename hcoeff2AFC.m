function  [varargout] = hcoeffRating(rawdata,mfr,nstims,rsp,fast, shiftwindows, bslshuff,shiftdur, varargin)

% old outputs kern,hcoeffs,hcoeffs2D
% 
% for outputs, include in order the kern, hcoeffs, hcoeffs_2d for each
% input, so if 2 inputs, 2 outputs in this order 
% varargout{1}=(kern1, hcoeff1, hcoeff2D1), varargout{2}= (kern2, hcoeff2, hcoef2D2)
%
%

% *********************************************************************************************
% *********************************************************************************************
% *******   This code was authored by Michael Hill (buffalohill@me.com), it may not     *******
% *******   be used as a whole or in parts without the written permission of Michael    *******
% *******   Hill and any use of this code (in whole or in parts) must reference the     *******
% *******   paper in which the h-coefficient method was first introduced (contact me    *******
% *******   to inquire about appropriate referencing formats for this code).             ******
% *********************************************************************************************
% *********************************************************************************************

% I guess peristim is the the activity in the responses to analyse, in a
% cell with each entry of the cell as a different response
% Doubt, must the perittls be already delimited by response time or just
% the raw one? 
% 

% nstims:
% nstims is the number of trials or stimuli that were run during the whole
% experiments
%
% peristim: in milliseconds, 
%
% rawdata in milliseconds
% rawdata is a vector containng timestamps for each spike that was recorded
% during yoru experiment.
%
% mfr in Hz
% mfr is the mean firing rate for this neuron across the whole
% experiment

% rsp are the start and end of the response time in milliseconds
%
% fast
% fast = 1 runs ranshuff in fast mode. fast = 0 runs randshuff in slow mode
%
% Output:
%
% Kern is the response sPSTH, so of all peristim responses, in Hz 
% 
% hcoeffs is a vector with the h coefficient for 1-positive deflection
% responses. 2- negative deflection responses
%
%hcoeffs2D contains the 3 parameters of h coefficient in this order (b,a,c)
% where b how much higher the response for the stimulus was than random, 
% a means how much wider the response was when compared to responses in
% random, and c is how high the random fluctuations where. Also for pos an
% neg deflactions. 
%


% bsl = [-1000,-200];
% rsp = [200,1000];
% nstims      =   1;
% mfr         =   3;
% fast        =   1;      	% Do the bootstrapped kernel method fast
%
% frombsl = round((bsl(1,1)/10) +500 + 1);
% tobsl = round((bsl(1,2)/10) +500);


dbstop if error




%% Generate a smooth PSTH from your response

out=1;
for k=1:length(varargin)
    eval('peristim = varargin{k};');
    
    
    for yuhu=1:shiftwindows(k)
        rsp1=rsp+(shiftdur*(yuhu-1));
        
        corrkern = ones(nstims,1)*NaN;
        
        for j = 1:nstims                                                            % Prep for the adaptive kernel analysis
            if j == 1
                forkern = peristim{1,1};
            else
                forkern = cat(2, forkern, peristim{j,1});
            end
            corrkern(j,1) = length(peristim{j,1});                                  % Preparation of the factor by which you need to multipy the ssvkern output, as it is normalized to the mean firing rate
        end
        
        [kern(2,:),kern(1,:)] = ssvkernel(forkern.*10^(-3),(-4.99:0.01:5));         % Smooth your PSTH
        
        kern(2,:) = kern(2,:) * (squeeze(nanmean(corrkern)));                       % Reverse the normaliztion applied by ssvkernel and transform to Hz.
        
        if isempty(kern(2,~isnan(kern(2,:))))                                       % Check for the presence of any firing rate whatsoever
            kern(2,:) = zeros(1,length(kern));
        end
        
        %% Find your local extrema
        
        fromrsp = round((rsp1(1,1)/10) +500 + 1);                                    % Find the index of the local max and local min of the response for further analysis
        torsp = round((rsp1(1,2)/10) +500);
        
        indmax = find(kern(2,fromrsp:torsp) == max(kern(2,fromrsp:torsp)),1,'first') + (fromrsp - 1);
        indmin = find(kern(2,fromrsp:torsp) == min(kern(2,fromrsp:torsp)),1,'first') + (fromrsp - 1);
        
        if isempty(indmax) || isempty(indmin)
            indmax = NaN;
            indmin = NaN;
        end
        
        %% Calculate the h-coefficient
        
        [hcoeff_max,hcoeff_min,hcoeff_max2D,hcoeff_min2D] = hcoeff_gen(bslshuff,kern(2,:),indmin,indmax,mfr);
        hcoeffs = [hcoeff_max,hcoeff_min];
        hcoeffs2D = [hcoeff_max2D;hcoeff_min2D];
        
        
        output{1}= kern;
        output{2}= hcoeffs;
        output{3}= hcoeffs2D;
        
        varargout{out}= output;
        out=out+1;

    end
end

end





function [hcoefficient_max,hcoefficient_min,hcoefficient_max2D,hcoefficient_min2D] = hcoeff_gen(bslshuff,kern,indmin,indmax,normit)
maxes = zeros(10001,1);                                                     % maxes & minses are two long vector of zeros, into which we can then enter the area measurements during the h-coefficient calculation further down
minses = maxes;
kern = (kern/normit)-1;                                                     % Normalize the response to the channel's baseline value and subtracts 1 so that the baseline value is = 0
if ~isnan(indmax) && kern(1,indmax) > 0                                     % This makes sure the maximum is above the mean firing rate of this channel, as else this analysis makes no sense
    firstdec = norm(floor(norm(kern(1,indmax))*10)/10);                     % Finds out what the next decimal point is between the local max and the channel baseline (in units of channel baseline firing rate, e.g. if the max firing rate is 1.568 times the mean firing rate of the channel, then your next decimal point towards that mean firing rate is 1.5)
    u = indmax(1); v = u;
    for b = 1:round(firstdec*10)+1
        while u > 1 && kern(1,u) > roundn(firstdec-((b-1)*0.1),-1)          % This travels along the curve towards the left (heading downwards from the local maximum) until it hits the next decimal multiple of the mean firing rate of the whole channel
            u = u-1;
        end
        while v < length(kern) && kern(1,v) > roundn(firstdec-((b-1)*0.1),-1) % This travels along the curve towards the right (heading downwards from the local maximum) until it hits the next decimal multiple of the mean firing rate of the whole channel
            v = v+1;
        end
        maxes(end-(round(firstdec*10)+1)+b,1) = trapz((kern(1,u:v))-roundn(firstdec-((b-1)*0.1),-1)).*0.1; % "-firstdec" makes sure that only the area between the next decimal multiple of the baseline firing rate and the peack that reaches just above this threshold is calculated. The factor 0.1 is there because the spacing of our x-axis is 0.1, see help trapz for more on this
    end
end
if ~isnan(indmin) && kern(1,indmin) < 0                                     % This makes sure the minimum is below the mean firing rate of this channel, as else this analysis makes no sense
    firstdec = norm(floor(norm(kern(1,indmin))*10)/10);                     % Finds out what the next decimal point is between your min and the channel baseline (in units of channel baseline firing rate, e.g. if the max firing rate is -0.568 times the mean firing rate of the channel, then your next decimal point towards that mean firing rate is 0.5)
    u = indmin(1); v = u;
    for b = 1:round(firstdec*10)+1
        while u > 1 && kern(1,u) < -roundn(firstdec-((b-1)*0.1),-1)         % This travels along the curve towards the left (heading upwards from the local minimum) until it hits the next decimal multiple of the mean firing rate of the whole channel
            u = u-1;
        end
        while v < length(kern) && kern(1,v) < -roundn(firstdec-((b-1)*0.1),-1) % This travels along the curve towards the right (heading upwards from the local minimum) until it hits the next decimal multiple of the mean firing rate of the whole channel
            v = v+1;
        end
        minses(end-(round(firstdec*10)+1)+b,1) = norm(trapz(kern(1,u:v)+roundn(firstdec-((b-1)*0.1),-1)).*0.1); % "+(firstdec-((b-1)*0.1))" makes sure that only the area between the next decimal multiple of the baseline firing rate and the peack that reaches just below this threshold is calculated. The factor 0.1 is there because the spacing of our x-axis is 0.1, see help trapz for more on this
    end
end
maxes = maxes(2:end,1) - maxes(1:end-1,:);                                  % From each row this subtracts the row right above it. Before doing so the area measurements represent the area above each decimal multiple of the baseline firing rate. After doing this, the area measurements are now corrected to reflect only the area inbetween this and the next following decimal. In other words, instead of reporting how large the area is above e.g. 1.7 * meanfr(channel), we now are reporting how large the area is inbetween 1.7 * meanfr(channel) and 1.8 * meanfr(channel), i.e. in each band
minses = minses(2:end,1) - minses(1:end-1,:);
hcoefficient_max = sum(maxes > bslshuff.maxes)/sum(bslshuff.maxes > 0);     % This now first counts in how many bands the response shows a bigger area (i.e. is wider) then the max in bands across the shuffstrapped samples. It then divides this number by the number of bands that showed any response at all in the shuffstrapped samples (i.e. the maximum reach of the highest baseline response). The resulting number means the following: If it is smaller then 1, then this means that the response was wider then the baseline responses in only some of the bands (although it may overall reach higher then any other, in which case hcoefficient < 1 would mean that at it's base this high response waas still narrower then the maxes from the baselines). If hcoefficient = 1 this means that the response was the same hight as the highest response in the baseline but wider in every band. If hcoefficient > 1 this means that the response was both higher then the baseline reseponses and wider in more bands.Any hcoefficient above zero can be reagegarded as significant dependent on your criteria. hcoefficient is higher then one should definitively be regarded as highly significant
hcoefficient_min = sum(minses > bslshuff.minses)/sum(bslshuff.minses > 0);
prov = (nnz(maxes)-nnz(bslshuff.maxes));                                    % number of bands that reac higher then in the shuffstrapped baseline
if prov > 0
    hcoefficient_max2D(1,1) = prov;
else
    hcoefficient_max2D(1,1) = 0;
end
prov = (nnz(maxes)-nnz(bslshuff.maxes));                                    % number of bands that reac higher then in the shuffstrapped baseline
if prov > 0
    hcoefficient_max2D(1,1) = prov;
else
    hcoefficient_max2D(1,1) = 0;
end
hcoefficient_max2D(1,2) = sum(maxes > bslshuff.maxes) - hcoefficient_max2D(1,1);  % number of bands that aren't higher but are wider then the in the shuffled baseline
hcoefficient_max2D(1,3) = sum(bslshuff.maxes > 0);                          % number of bands in the shuffled baseline
prov = (nnz(minses)-nnz(bslshuff.minses));                                  % number of bands that reach higher then in the shuffled baseline
if prov > 0
    hcoefficient_min2D(1,1) = prov;
else
    hcoefficient_min2D(1,1) = 0;
end
hcoefficient_min2D(1,2) = sum(minses > bslshuff.minses) - hcoefficient_min2D(1,1);  % number of bands that aren't higher but are wider then the in the shuffled baseline
hcoefficient_min2D(1,3) = sum(bslshuff.minses > 0);                         % number of bands in the shuffled baseline
end


