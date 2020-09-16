function [ imageout ] = ColourSuperpixels( labels, RGBs )
%COLOURSUPERPIXELS Takes as input a matrix of labels and corresponding
%image and assigns each superpixel its average colour, then returns an RGB
%matrix

reds = uint8(labels*0);
blues = uint8(labels*0);
greens = uint8(labels*0);

[h, w, ~] = size(labels);

for y = 1:h
    for x = 1:w
        label = labels(y,x);
        reds(y,x) = RGBs(label+1,1);
        greens(y,x) = RGBs(label+1,2);
        blues(y,x) = RGBs(label+1,3);
    end
end


imageout = cat(3,reds,greens);
imageout = cat(3,imageout,blues);

end

