function [satcue, ranking, events, onset_times, DirName] = load_behavioural_data (DirName)

% this function load the data from a specific patient
% 

dbstop if error
if nargin ==0
    [path]= uigetdir('/media/Projects/Alex/Satiation Analysis', 'Select a patient folder');
    [~,DirName]=fileparts(path);
else
    path=DirName;
end
cd (path)



%load behavioural results 
pathfolders=dir;
for i=1:length(pathfolders)
    if regexp(pathfolders(i).name, 'Files')
       cd (pathfolders(i).name)
       filesfolder=dir;
       for j=1:length(filesfolder)
         if regexp(filesfolder(j).name, 'new_sat_output')
             load (filesfolder(j).name)
         end
         if regexp(filesfolder(j).name, 'alexevents')
             load (filesfolder(j).name)
         end
       end
    end
    
end
cd ../..
end


