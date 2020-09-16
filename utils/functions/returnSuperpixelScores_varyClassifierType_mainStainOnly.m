function [ pixelScores ] = returnSuperpixelScores_varyClassifierType_mainStainOnly(classifierPath, labels, modelType,features)
%classifySuperpixelsInTile classifies the superpixels in a given tile by
%generating a feature vector and applying the classifiers in SVMModels

load(classifierPath);
% This loads in "classes" and "SVMModels" or similar
if strcmp(modelType,'SVM') % Support Vector Machine
    models = SVMModels;
elseif strcmp(modelType,'NB') % Naive Bayes
    models = NBModels;
elseif strcmp(modelType,'KNN') % K nearest neighbours
    models = KNNModels;
end
Features = features';

numRequiredFeaturesForClassifier = length(models{2,1}.PredictorNames);
if size(Features,2) > numRequiredFeaturesForClassifier
   Features = Features(:,1:numRequiredFeaturesForClassifier); 
end

% Generate SVMs
% predictedClassification = zeros(length(Features),1);

N = size(Features,1);
Scores = zeros(N,1);

% for j = 2%1:numel(classes)
    [~,score] = predict(models{2},Features);
    [~,w] = size(score);
    if w == 1
        Scores(:) = score(:,1);%2); % Second column contains positive-class scores
    else
        Scores(:) = score(:,2);%2); % Second column contains positive-class scores
    end
% end
% [~,maxScoreIndexr] = max(Scores,[],2);
% BinaryScores = Scores > 0;

% Colour in corresponding labels in prediction - colour by cell type
[h, w, ~] = size(labels);
pixelScores = zeros(h,w);

for i = 1:h
    label = labels(i,:);
    pixelScores(i,:) = Scores(label);
end

% for i = 1:h
%     for k = 1:3
%         label = labels(i,:);
%         pixelScores(i,:,k) = Scores(label,k);
%     end
% end


end

