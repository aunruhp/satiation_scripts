function [correlations_rating_ranking] = correlation_between_rating_ranking (ranking)

%{ 
This function sorts and compares the ordering of preference during 
satiation paradigm according to 3 parameters:1-added value given during rating 
task (from -300 to 300); number of choices of the stimulus during 2AFC task
(from 0 to 19), and added value given in 2AFC task (from 0 to 1900).
In order to check for the consistency of this measures. 

Input: rating is a variable saved during Satiation paradigm, saved in 
runsat_output_*.mat
It contains the information about which value was given to each stimulus
during the rating/evaluation task and during the 2AFC task.


%}

    PRE_sorted_by_rating= sortrows(ranking{2}, 3);               %sorted by rating
    PRE_sorted_by_ranking= sortrows(ranking{3}, 3);     %sorted by ranking during 2AFC
    PRE_sorted_by_rating2AFC= sortrows(ranking{3}, 4);  %sorted by rating during 2AFC

    POST_sorted_by_rating= sortrows(ranking{5}, 3);     %sorted by rating
    POST_sorted_by_ranking= sortrows(ranking{6}, 3);    %sorted by ranking during 2AFC
    POST_sorted_by_rating2AFC= sortrows(ranking{6}, 4); %sorted by rating during 2AFC

    combinations=nchoosek([{'by_rating'}, {'by_ranking'}, {'by_rating2AFC'}], 2);

    correlations_rating_ranking=cell(1,8);
    correlations_rating_ranking{1}='Before meal, 2st rating vs ranking, 3nd rating vs rating2AFC, 4rd ranking vs rating2AFC';
    correlations_rating_ranking{5}='After meal,  6st rating vs ranking, 7nd rating vs rating2AFC, 8rd ranking vs rating2AFC';

    for i=1:3
        eval(['correlations_rating_ranking{i+1}=corr(PRE_sorted_' combinations{i,1} '(:,1), PRE_sorted_' combinations{i,2} '(:,1));']);   
        eval(['correlations_rating_ranking{i+5}=corr(POST_sorted_' combinations{i,1} '(:,1), POST_sorted_' combinations{i,2} '(:,1));']); 
    end




end