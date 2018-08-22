ObjDir = 'C:/Users/mxxh/Desktop/singal project/source code/src/samples/';
beta_vector = zeros(2,2);
imname = dir('C:/Users/mxxh/Desktop/singal project/source code/src/samples/*.tif');
im_num = length(imname);
for i = 1:1:im_num  
    bgFile = [ObjDir,imname(i).name];
    beta = extract_feature(bgFile);
    beta_vector(i, :) = beta;
end  