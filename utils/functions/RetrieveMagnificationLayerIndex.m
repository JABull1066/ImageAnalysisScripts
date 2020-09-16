function [ i ] = RetrieveMagnificationLayerIndex( path, magnification )
%RetrieveMagnificationLayerIndex From an image (NDPI or SVS), retrieve the specified
%magnification
disp(path)

if strcmp(path(end-3:end),'ndpi')
    info = imfinfo(path);

    mag = 0;
    for i = 1:length(info)
        % Loop tag IDs until we find ID 65421
        for j = 1:length(info(i).UnknownTags)
            ID = info(i).UnknownTags(j).ID;
            if ID == 65421
                mag = info(i).UnknownTags(j).Value;
                break
            end
        end
        if mag == magnification
            break
        end
    end

    if mag ~= magnification
        message = ['Magnification ',num2str(magnification),' not found'];
        error(message)
    end
elseif strcmp(path(end-2:end),'svs')
    info = imfinfo(path);

    imageDetails = split(info(1).ImageDescription,'|');
    maxWidth = info(1).Width;
    for j = 1:length(imageDetails)
        if strcmp(imageDetails{j}(1:6),'AppMag')
            temp = split(imageDetails{j},' ');
            maxMag = str2num(temp{3});
        end
    end
    
    mag = 0;
    for i = 1:length(info)
        mag = round(info(i).Width * maxMag / maxWidth,3);
        if mag == magnification
            break
        end
    end
    
    if mag ~= magnification
        message = ['Magnification ',num2str(magnification),' not found'];
        error(message)
    end
elseif strcmp(path(end-2:end),'scn')
    % Designed for use on the Manchester HN images. This is a bit of a
    % bodge.
    % ASsume images scanned at 20x mag (because they were)
    info = imfinfo(path);
    
    % Find second occurence of r="0" - first occurence describes high res
    str = info(1).ImageDescription;
    startInds = strfind(str,'r="0"');
    endInds = strfind(str(startInds(2):end),' />');
    ifd = str(startInds(2)+6:startInds(2)+endInds(1)-2);
    temp = split(ifd,'"');
    layer = str2num(temp{2})+1; %As Matlab indexes from 1, not 0
    % OK, that gives us the layer index of the highest res actual picture.
       
    maxMag = 20;
    mag = 0;
    for i = layer:length(info)
        mag = round(info(i).XResolution * maxMag / info(layer).XResolution,3);
        if mag == magnification
            break
        end
    end
    
    if mag ~= magnification
        message = ['Magnification ',num2str(magnification),' not found'];
        error(message)
    end
else
    error('File type not supported for RetrieveMagnificationLayerIndex - is it NDPI or SVS?')    
end

end

