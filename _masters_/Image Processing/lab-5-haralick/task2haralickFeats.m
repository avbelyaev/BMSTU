
files = dir('*.png');
IMAGE_NUM = numel(files);
HARALICK_FEAT_NUM = 14;

haralicks = zeros(IMAGE_NUM, HARALICK_FEAT_NUM);
i = 1;

for file = files'
    S=imread(file.name);
    S=rgb2gray(S);
    I= imresize (S, [350 350]);
    glcm=graycomatrix(I,'offset',[-1 1],'NumLevel', 8,'Symmetric',true);
    xFeatures = 1:HARALICK_FEAT_NUM;
    
    j = 1;
    feats = haralickTextureFeatures(glcm, xFeatures);
    while j <= HARALICK_FEAT_NUM
        haralicks(i, j) = feats(j);
        j = j + 1;
    end
    i = i + 1;
end
NUM_CLUSTERS = 4;
%% K-means
% clusterized = kmeans(haralicks, NUM_CLUSTERS);

%% Hierachial
Y = pdist(haralicks);
Z = linkage(Y);
dendrogram(Z);

clusterized = cluster(Z, 'maxclust', NUM_CLUSTERS);

disp(clusterized);
