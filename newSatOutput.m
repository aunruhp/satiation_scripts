function newSatOutput

% This function loads the output from satiation for the chosen patient
% corrects its sign (a bug was found in the paradigm code) and reckons a
% new variable, namely, a new lukert 2AFC output, from -1900 to +1900,
% which contains not only the responses when an object is chosen, but also
% the negative ones and adds them all. This is saved in a new .mat file
% called new_runsat_output_p0??.mat which contains corrected satcue and
% corrected and expanded ranking with the new variable. 


dbstop if error
[path]= uigetdir('/media/Projects/Alex/Satiation Analysis', 'Select a patient folder');
cd (path)
%load behavioural results
pathfolders=dir;
for i=1:length(pathfolders)
    if regexp(pathfolders(i).name, 'Files')
        cd (pathfolders(i).name)
        filesfolder=dir;
        for j=1:length(filesfolder)
            if regexp(filesfolder(j).name, 'runsat_output')
                load (filesfolder(j).name)
                break
            end
        end
    end
    
end


%cd (path)

%correct the errors in the sign of rating
satcue{1,2}(:,14) = mat2cell(cellfun(@(x) -x,satcue{1,2}(:,14)),ones(1,size(satcue{1,2},1)),1);
satcue{1,5}(:,14) = mat2cell(cellfun(@(x) -x,satcue{1,5}(:,14)),ones(1,size(satcue{1,5},1)),1);
% Fix ranking (inverted Likert in rating trials):
ranking{1,2}(:,3) = 3-ranking{1,2}(:,3);
ranking{1,5}(:,3) = 3-ranking{1,5}(:,3);
ranking{1,2}(:,4) = -ranking{1,2}(:,4);
ranking{1,5}(:,4) = -ranking{1,5}(:,4);

newVal2AFC = nan(20,2);

for i=1:20

    
    ind3l= find(cell2mat(satcue{3}(:,2))==i);
    ind3r= find(cell2mat(satcue{3}(:,3))==i);
    ind4l=find(cell2mat(satcue{6}(:,2))==i);
    ind4r=find(cell2mat(satcue{6}(:,3))==i);
    

   
    newVal2AFC(i,1)= sum([-(cell2mat(satcue{3}(ind3l,14)));(cell2mat(satcue{3}(ind3r,14)))]);
    newVal2AFC(i,2)= sum([-(cell2mat(satcue{6}(ind4l,14)));(cell2mat(satcue{6}(ind4r,14)))]);
    
end


ranking{3}(:,5)= newVal2AFC(:,1);
ranking{6}(:,5)= newVal2AFC(:,2);

save (['new_' filesfolder(j).name(4:end)], 'satcue',  'ranking')

end