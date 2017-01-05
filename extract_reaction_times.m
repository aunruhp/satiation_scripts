function [RT] = extract_reaction_times(events)

% This function extracts the real reaction times of satiation, 
% i.e., the time between the likert is shown to the first 
% press buttom
% 
% Input: events of satiation. Events contains the timestamps of each event
% in the first column and the event_ID in the second
% Output: reaction times

ind1=find(events(:,2)==50);
 

if all(events(ind1-1,2)==44) &&  length(ind1)==120
    RT{2}=events(ind1(1:60)+1,1) - events(ind1(1:60),1); 
    RT{5}=events(ind1(61:end)+1,1) - events(ind1(61:end),1); 

else 
    error ('Event 44 should be followed by event 50 in rating trials and 120 trials, some error here!')
end


ind3=find(events(:,2)==51);

if all(events(ind3-1,2) == 43) && length(ind3)== 380
    RT{3}=events(ind3(1:190)+1,1) - events(ind3(1:190),1); 
    RT{6}=events(ind3(191:end)+1,1) - events(ind3(191:end),1); 

else 
    error ('Event 43 should be followed by event 51 in 2AFC trials and 380 trials, some error here!')
end

end

