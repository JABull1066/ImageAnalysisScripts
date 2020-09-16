function [micronsPerPixel, pixelsPerMm] = GetMicronsPerPixel(pathToImage)

temp = split(pathToImage,'.');
filetype = temp{end};
% if .SVS
micronsPerPixel = 0;
if strcmp(filetype,'svs')
    info = imfinfo(pathToImage);
    details = split(info(1).ImageDescription,'|');
    for i = 1:length(details)
        if strcmp(details{i}(1:3),'MPP')
            temp = split(details{i},' ');
            micronsPerPixel = str2num(temp{3});
        end
    end
    if micronsPerPixel == 0
       error('Could not identify microns per pixel for .svs image') 
    end
elseif strcmp(filetype,'ndpi')
    NDPI_info = imfinfo(pathToImage);
    if ~strcmp(NDPI_info(1).ResolutionUnit,'Centimeter')
       error('NDPI resolution should be "Centimeter"')
    end
    micronsPerPixel=10000/NDPI_info(1).XResolution;
    if micronsPerPixel == 0
       error('Could not identify microns per pixel for .ndpi image') 
    end
elseif strcmp(filetype,'scn')
    % Designed for use on the Manchester HN images. This is a bit of a
    % bodge.
    info = imfinfo(pathToImage);
    
    % Find second occurence of r="0" - first occurence describes high res
    str = info(1).ImageDescription;
    startInds = strfind(str,'r="0"');
    endInds = strfind(str(startInds(2):end),' />');
    ifd = str(startInds(2)+6:startInds(2)+endInds(1)-2);
    temp = split(ifd,'"');
    layer = str2num(temp{2})+1; %As Matlab indexes from 1, not 0
    % OK, that gives us the layer index of the highest res actual picture.
    micronsPerPixel=10000/info(layer).XResolution;
else
    error('File type not recognised - is it .svs, .scn, or .ndpi?')
end

pixelsPerMm = (1e3)/micronsPerPixel;

end

