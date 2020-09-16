function [ bw3, centroids ] = IdentifyCellsViaWatershed_DependentOnSize( mask, image, MPP,cellDiameterRangeInMicrons)
%IDENTIFYCELLSVIAWATERSHEDAUGMENTED Take in a layer of a mask (e.g., binary
%identified macrophage classification) and identify cell locations
%
% See https://blogs.mathworks.com/steve/2013/11/19/watershed-transform-question-from-tech-support/
bw = mask;

% Default
% excludeArea = 1000;
% assume cell area is 20 mum x 20 mum = 400 mum
% Exclude anything smaller than 10 mum x 10 mum
lowerSizeThreshold = cellDiameterRangeInMicrons(1);
upperSizeThreshold = cellDiameterRangeInMicrons(2);

% Cells smaller than this will be removed
excludeArea = round(lowerSizeThreshold * MPP * lowerSizeThreshold * MPP);
bw2 = bwareaopen(bw, excludeArea);
% bw2 = ~bwareaopen(~bw, 10);
% imshow(bw2)

D = -bwdist(~bw);
% imshow(D,[])

% Ld = watershed(D);
% imshow(label2rgb(Ld))

% bw2 = bw;
% bw2(Ld == 0) = 0;
% imshow(bw2)

mask2 = imextendedmin(D,2);
% imshowpair(bw,mask2,'blend')

D2 = imimposemin(D,mask2);
Ld2 = watershed(D2);
bw3 = bw;
bw3(Ld2 == 0) = 0;
% Also remove cells excluded by the small area removal in bw2
bw3(bw2 == 0) = 0;
% imshow(bw3)

% Cells larger than this will be split into two
divideArea = round(upperSizeThreshold * MPP * upperSizeThreshold * MPP);


% For any labels which are too large, we will split the corresponding pixel
% in two. Find longest axis and divide perpendicularly (assuming by this
% point all connected components are roughly elliptical
allRegionsSorted = false;
counter = 0;
lastTime = intmax;
while ~allRegionsSorted && (counter < 100)
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
%                 disp(k)
            end
            % Want all pixels of bw3
            % ...laying along minor axis
            % ....inside L == i
            %        point1 =
            % Minor axes
            xbar = s(k).Centroid(1);
            ybar = s(k).Centroid(2);
            
            b = 2*s(k).MinorAxisLength/2; % extra factor of 1.5 to ensure "cell" is fully contained within ellipse
            
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
    debug = false;
%     if debug
%         imshow(bw3);
    if length(badsizes) == 1
        disp(['Size: ',num2str(badsizes(1))])
        if lastTime == badsizes(1)
            disp('Unable to reduce final region size further - continuing')
            allRegionsSorted = true;
        end
        lastTime = badsizes(1);
    end
%     end
end

% Remove small areas again
excludeArea = round(lowerSizeThreshold * MPP * lowerSizeThreshold * MPP);
bw3 = bwareaopen(bw3, excludeArea);

%
% % For visualising
% imshow(bw3)
% hold on
%
% phi = linspace(0,2*pi,50);
% cosphi = cos(phi);
% sinphi = sin(phi);
%
% for k = 1:length(s)
%     xbar = s(k).Centroid(1);
%     ybar = s(k).Centroid(2);
%
%     a = s(k).MajorAxisLength/2;
%     b = s(k).MinorAxisLength/2;
%
%     theta = pi*s(k).Orientation/180;
%     R = [ cos(theta)   sin(theta)
%          -sin(theta)   cos(theta)];
%
%     xy = [a*cosphi; b*sinphi];
%     xy = R*xy;
%
%     x = xy(1,:) + xbar;
%     y = xy(2,:) + ybar;
%
%     plot(x,y,'r','LineWidth',2);
%
%     % Minor axes
%     angle = theta;
%     plot([xbar-b*sin(angle),xbar+b*sin(angle)],[ybar-b*cos(angle),ybar+b*cos(angle)],'g')
%
% end
% hold off
%



s = regionprops(bw3,'centroid');
centroids = cat(1,s.Centroid);

end

