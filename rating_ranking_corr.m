function rating_ranking_corr

% this function analyses different aspects of the relationship of
%rating and ranking


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
         end
       end
    end
      
end




ranking{3}(:,4)
ranking{3}(:,3)
ranking{2}(:,4)



end