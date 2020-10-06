function [ bw3, centroids ] = IdentifyCellsViaWatershed_DependentOnSize( mask, image, MPP,cellDiameterRangeInMicrons)
%IDENTIFYCELLSVIAWATERSHEDAUGMENTED Take in a layer of a mask (e.g., binary
%identified macrophage classification) and identify cell locations
%
% See https://blogs.mathworks.com/steve/2013/11/19/watershed-transform-question-from-tech-support/
bw = mask;

lowerSizeThreshold = cellDiameterRangeInMicrons(1);
upperSizeThreshold = cellDiameterRangeInMicrons(2);

% Cells smaller than this will be removed
excludeArea = round(lowerSizeThreshold * MPP * lowerSizeThreshold * MPP);
bw2 = bwareaopen(bw, excludeArea);

D = -bwdist(~bw);

mask2 = imextendedmin(D,2);

D2 = imimposemin(D,mask2);
Ld2 = watershed(D2);
bw3 = bw;
bw3(Ld2 == 0) = 0;
% Also remove cells excluded by the small area removal in bw2
bw3(bw2 == 0) = 0;

% Cells larger than this will be split into two
divideArea = round(upperSizeThreshold * MPP * upperSizeThreshold * MPP);


% For any labels which are too large, we will split the corresponding pixel
% in two. Find longest axis of surrounding ellipse and divide perpendicularly
allRegionsSorted = false;
counter = 0;
lastTime = intmax;
while ~allRegionsSorted && (counter < 100) % Allow up to 100 iterations
    counter = counter + 1;
    disp(['Iteration ',num2str(counter),' of shrinking regions'])
    badRegions = 0;
    [B,L,N] = bwboundaries(bw3);
    s = regionprops(L, 'Area','Orientation', 'MajorAxisLength','MinorAxisLength', 'Eccentricity', 'Centroid');
    allRegionsSorted = true;
    badsizes = [];
    for k = 1:length(s)
        if s(k).Area >  divideArea
            if 0.5*s(k).Area > divideArea
                allRegionsSorted = false;
                badRegions = badRegions + 1;
                badsizes(badRegions) = 0.5*s(k).Area;
            end
            % Want all pixels of bw3
            % ...laying along minor axis
            % ....inside L == i
            % Minor axes
            xbar = s(k).Centroid(1);
            ybar = s(k).Centroid(2);
            
            b = 2*s(k).MinorAxisLength/2; % extra factor of 2 to ensure "cell" is fully contained within ellipse
            
            theta = pi*s(k).Orientation/180;
            angle = theta;
            x1 = xbar-b*sin(angle);
            x2 = xbar+b*sin(angle);
            y1 = ybar-b*cos(angle);
            y2 = ybar+b*cos(angle);
            %         plot([x1,x2],[y1,y2],'g')
            
            % Create a line from point 1 to point 2
            numSamples = ceil(sqrt((x2-x1)^2+(y2-y1)^2) / 0.4);
            x = linspace(x1, x2, numSamples);
            y = linspace(y1, y2, numSamples);
            xy = round([x',y']);
            dxy = abs(diff(xy, 1));
            duplicateRows = [0; sum(dxy, 2) == 0];
            % Now for the answer:
            finalxy = xy(~duplicateRows,:);
            finalx = finalxy(:, 1);
            finaly = finalxy(:, 2);
            for i = 1:length(finalx)
                try
                    if L(finaly(i),finalx(i)) == k
                        bw3(finaly(i),finalx(i)-1) = 0;% "blur" line slightly to ensure no 4-connectivity remaining
                        bw3(finaly(i)-1,finalx(i)) = 0;% "blur" line slightly to ensure no 4-connectivity remaining
                        bw3(finaly(i),finalx(i)) = 0;
                        bw3(finaly(i)+1,finalx(i)) = 0; % "blur" line slightly to ensure no 4-connectivity remaining
                        bw3(finaly(i),finalx(i)+1) = 0; % "blur" line slightly to ensure no 4-connectivity remaining
                        % More blurring
                        bw3(finaly(i)+1,finalx(i)+1) = 0;
                        bw3(finaly(i)-1,finalx(i)+1) = 0;
                        bw3(finaly(i)-1,finalx(i)-1) = 0;
                        bw3(finaly(i)+1,finalx(i)-1) = 0;
                    end
                catch
                    continue
                end
            end
        end
    end
    disp([num2str(badRegions),' large regions remaining'])

    if length(badsizes) == 1
        disp(['Size: ',num2str(badsizes(1))])
        if lastTime == badsizes(1)
            disp('Unable to reduce final region size further - continuing')
            allRegionsSorted = true;
        end
        lastTime = badsizes(1);
    end
end

% Remove small areas again
excludeArea = round(lowerSizeThreshold * MPP * lowerSizeThreshold * MPP);
bw3 = bwareaopen(bw3, excludeArea);

s = regionprops(bw3,'centroid');
centroids = cat(1,s.Centroid);

end

