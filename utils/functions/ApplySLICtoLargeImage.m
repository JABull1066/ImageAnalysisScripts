function [superpixelImage_master, labels_master, LABintensities_master, LABvariances_master, features_master] = ApplySLICtoLargeImage( originalImage,desiredSLICImagesize, desiredSuperpixelSize,superpixelColourSpaceWeightingParameter)
%ApplySLICtoLargeImage take a single image and superpixellate it
% If nothing else, superpixelColourSpaceWeightingParameter should equal 20
% by default 

%desiredSLICImagesize = 1000; % We want to SLIC only 1000 pixels at a time
overlapSize = 0;%25; % make the tiles overlap



% superpixelImage_master = uint8(nan(size(originalImage)));
% labels_master = int32(nan(size(originalImage(:,:,1))));
superpixelImage_master = zeros(size(originalImage),'uint8');
labels_master = zeros(size(originalImage(:,:,1)),'int32');

[h,w,~] = size(originalImage);

hStarts = [1:desiredSLICImagesize:h];
hEnds = hStarts + desiredSLICImagesize - 1;
hEnds(end) = h;
wStarts = [1:desiredSLICImagesize:w];
wEnds = wStarts + desiredSLICImagesize - 1;
wEnds(end) = w;

% Number of times we must apply SLIC algorithm
%numberOfPartitions = length(hStarts)*length(wStarts);
%desiredSuperpixelSize = 15;%30;
maxLabel = 0; % counter for label uniqueness
LABintensities_master = [];
LABvariances_master = [];
features_master = [];

for slicNumberH = 1:length(hStarts)
    for slicNumberW = 1:length(wStarts)
        %disp([slicNumberH, slicNumberW])
        % Now apply SLIC algorithm to each section in turn
        hStart = hStarts(slicNumberH) - overlapSize;
        wStart = wStarts(slicNumberW) - overlapSize;
        hEnd = hEnds(slicNumberH) + overlapSize;
        wEnd = wEnds(slicNumberW) + overlapSize;
        % Deal with cases where we on't want an overlap - i.e., edges
        if slicNumberH == 1
            hStart = hStarts(slicNumberH);
        end
        if slicNumberW == 1
            wStart = wStarts(slicNumberW);
        end
        if slicNumberH == length(hStarts)
            hEnd = hEnds(slicNumberH);
        end
        if slicNumberW == length(wStarts)
            wEnd = wEnds(slicNumberW);
        end
        
        % Now we SLIC the image tile that we've found
        image = originalImage(hStart:hEnd,wStart:wEnd,:);
        [hWid,wWid,~] = size(image);
        NumPixels = round(hWid*wWid/desiredSuperpixelSize);
        
        %[labels, numlabels, LABintensities, ~, seedIndex, LABvariances] = slicmexExtendedFeaturesJB(image,NumPixels,20);%numlabels is the same as number of superpixels
        %[labels, numlabels, LABintensities, ~, seedIndex, LABvariances, neighbouringLabelsStore, neighbouringLabelsCounter] = slicmexExtendedFeaturesAndReturnNeighboursJB(image,NumPixels,20);%numlabels is the same as number of superpixels
        %[labels, numlabels, LABintensities, ~, seedIndex, LABvariances, neighbouringLabelsStore, neighbouringLabelsCounter, twoStepNeighbouringLabelsStore, twoStepNeighbouringLabelsCounter] = slicmexExtendedFeaturesAndReturnNeighboursJB_backupTemp(image,NumPixels,20);%numlabels is the same as number of superpixels
%         [labels, numlabels, LABintensities, ~, seedIndex, LABvariances, features] = slicmexReturnAllFeatures(image,NumPixels,20);%numlabels is the same as number of superpixels

        [labels, numlabels, LABintensities, ~, seedIndex, LABvariances, features] = slicReturnExtendedFeatures(image,NumPixels,superpixelColourSpaceWeightingParameter);%numlabels is the same as number of superpixels
        
%         % OK, this one returns 26 features. We only want 18 for this
%         % classifier.
%         features = features(1:18,:);
%         
        % Colour superpixels
        LABintensitiesPrime = single(LABintensities');
        RGBintensities = lab2rgb(LABintensitiesPrime,'OutputType','uint8');
        [superpixelImage] = ColourSuperpixels( labels, RGBintensities );
        superpixelImage_master(hStart:hEnd,wStart:wEnd,:) = superpixelImage(:,:,:);
        labels_master(hStart:hEnd,wStart:wEnd) = labels + maxLabel;
%        neighbouringLabelsStore = neighbouringLabelsStore + maxLabel;
        maxLabel = maxLabel + numlabels;
        
        LABintensities_master = cat(2,LABintensities_master,LABintensities);
        LABvariances_master = cat(2,LABvariances_master,LABvariances);
        features_master = cat(2,features_master,features);
    end
end
      
% For simplicity, we add one to every label so that our indexing starts
% from 1 instead of 0.
labels_master = labels_master + 1;
end

