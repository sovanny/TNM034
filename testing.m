%% for testing stuff
clc
clear
close all
images = dir('./Images_Training/*.jpg');

%%
result = [];
for img = images'
    img
    %result = [result,tnm034(strcat(img.folder,'\',img.name))];
    result = [result,tnm034(strcat(img.folder,'/',img.name))];
    
end

%%
tnm034('./Images_Training/Le_1_Example.jpg')
