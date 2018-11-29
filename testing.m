%% for testing stuff
images = dir('./Images_Training/*.jpg');
%%
result = [];
for img = images'
    
    %result = [result,tnm034(strcat(img.folder,'\',img.name))];
    result = [result,tnm034(strcat(img.folder,'/',img.name))];
    
end