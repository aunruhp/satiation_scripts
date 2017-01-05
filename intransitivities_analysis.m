function [pval, intransitivity_ratio, diffValue]=intransitivities_analysis(DirName)

% This fuction analyses the existence of significant intransitivity in the 
% decisions made by a patient. For this, it takes the rating value of each
% of the products, compares the value of both in 2AFC, and checks if the
% direction of the decision is intransitive in a binomial test with 0.05 of
% thershold. 
%
% OUTPUT:
% pval is the result of the binomial test for presence of intransitivities
% intransitivity_ratio is the ratio of intransitivities in this patient
% diff_intran is the absolute value difference when the intransitivity
% occured




if nargin ==0
    [satcue, ranking, events] = load_behavioural_data;
elseif nargin>0
     [satcue, ranking, events] = load_behavioural_data(DirName);
end



%% Find the identities of left and right images during 2AFC

k=1;
for jj=[3,6]
    
ident_left= cell2mat(satcue{jj}(:,2));
ident_right= cell2mat(satcue{jj}(:,3));

%% Create looping variables and 


j=1;
for i= ident_left'
    rank1(j,1)= ranking{2}(i,4); %for each left image, check ID and include the ranking obtained during evaluation
    j=j+1;
end


j=1;
for i= ident_right'
    rank2(j,1)= ranking{2}(i,4);
    j=j+1;
end


rank_difference=rank1>rank2;
value=cell2mat(satcue{jj}(:,14));
diff_intran=[];
intrans=0;
h=1;
for i=1:length(rank_difference)
    if rank_difference(i)==1                                %if left higher value than right, value should be negative, which means left chosen
        if value(i)>0
%             disp(['intransitivity found in ' num2str(i)]);
            intrans=intrans +1;
            diff_intran(h) = abs(rank1(i)-rank2(i));
            h=h+1;
        end
    end
    if rank_difference(i)==0                                %if RIGHT higher value than LEFT, value should be POSITIVE, which means RIGHT chosen
        if value(i)<0
%             disp(['intransitivity found in ' num2str(i)]);
            intrans=intrans +1;
            diff_intran(h) = abs(rank1(i)-rank2(i));
            h=h+1;
        end
    end
end

intransitivity_ratio{k}=intrans/length(ident_left);


pval{k} = myBinomTest(intrans,length(ident_left),0.05,'one');

diffValue{k}= diff_intran;

k=k+1;


end



% %biomial test
% pval= sum(binopdf(intrans:length(ident_left), length(ident_left), 0.05));
% 
% 
% %check this!
% if intrans <= 190/2
%     p_bino_l = sum(binopdf([1:intrans], 190, 0.5));
%     if p_bino_l< 0.05
% %     disp('This patient has a significant intransitivity')
%     end
% elseif intrans > 190/2
%     p_bino_r = sum(binopdf([intrans:190], 190, 0.5));
%     if p_bino_r<0.05
% %     disp('This patient has a significant intransitivity')
%     end
% end
