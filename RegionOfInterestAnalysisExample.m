% This is a demonstration script for analysing a region of interest from an
% IHC slide. This code uses SLIC superpixellation (see, e.g.
% https://www.epfl.ch/labs/ivrl/research/slic-superpixels/) to divide the
% image into superpixels, and then uses a support vector machine (SVM)
% classifier to identify cells. Finally, approximate cell centres are
% identified via a modified watershedding alogirthm.

pathToImage = 'ExampleTif_CD68.tif';
pathToClassifier = 'Classifiers/ExampleClassifierCD68.mat';

micronsPerPixel = 0.882;

cellDiameterRangeInMicrons = [8,21];

img = imread(pathToImage);
[h,w,d] = size(img);

% Parameters for SLIC analysis
desiredSLICImagesize = 1000; % Size of square region on which to apply SLIC - useful for large images. Region will be analysed in squares of side length desiredSLICImagesize
desiredSuperpixelSize = 20; % Desired superpixel size in pixels
superpixelColourSpaceWeightingParameter = 20; % Parameter

[superpixels, labels, LABint, LABvar, features] = ApplySLICtoLargeImage( img,desiredSLICImagesize, desiredSuperpixelSize,superpixelColourSpaceWeightingParameter);
[ scores ] = returnSuperpixelScores_varyClassifierType_mainStainOnly(pathToClassifier, labels, 'SVM', features);
immuneStainThreshold = 0;
mask = scores(:,:)>immuneStainThreshold;

[ abw2, immuneCentroids ] = IdentifyCellsViaWatershed_DependentOnSize( mask, img, micronsPerPixel,cellDiameterRangeInMicrons);

imshow(img)
hold on
plot(immuneCentroids(:,1),immuneCentroids(:,2),'g.','LineWidth',2,'MarkerSize',10)
hold off