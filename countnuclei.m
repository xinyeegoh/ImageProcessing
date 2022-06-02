function result=countnuclei(imageinput)

%%
% function INPUT argument: imageinput, which is the root cell image
% function OUTPUT argument: result, which is a structure type storing
%                                   neccesary variables 

%%
%PRE-STEP: make a copy of the original input image before processing
result.img=imageinput; 

%%
%STEP 1: remove RED and BLUE components from the image
%why?
%   because we are only detecting the nuclei which are the green components in the root cell image
    
%image is RGB type, so (:,:,:1) represents RED (:,:,3) represents BLUE
result.img(:,:,1)=0; %REMOVE RED
result.img(:,:,3)=0; %REMOVE BLUE

figure;
subplot(2,2,1);
imshow(result.img); title('GREEN extracted from input image');


%%
%STEP 2: convert image from RGB to GRAY-SCALE
%since we have previously extracted only the GREEN COMPONENTS 
%gray intensity shall take effective on the GREEN COMPONENTS only 
%and the rest shall be black or 0 intensity
result.img_g = rgb2gray(result.img);
[v,h]=size(result.img_g);

subplot(2,2,2);
imshow(result.img_g); title('GRAY-SCALE of image');


%%
%STEP 3: adjust image to increase contrast
%why?
%   so that the nuclei can stand out more from the dark surrounding
%   and we can observe the nucleii clearer
    
result.img_c = imadjust(result.img_g);

subplot(2,2,3);
imshow(result.img_c); title('increase image contrast');


%%
%STEP 4: remove noise from the image through smoothing filter

%using median filter because it is good at retaining the image contrast
result.img_c=medfilt2(result.img_c,[3 3],'zeros');

subplot(2,2,4);
imshow(result.img_c);title('smooth image for noise reduction');


%%
%STEP 5a: edge detection of the smoothened image
%       to detect the boundaries of the nuclei

%using CANNY edge detection, why?
%   it is good at not getting fooled by noises as "edges"
%   it can detect weak edges which are the less "green" nuclei in this case

%detect using canny edge detection 
result.edge=edge(result.img_c,'canny',[],2);

figure; sgtitle('Edge Detection + Processing');
subplot(3,1,1);
imshow(result.edge);title('image canny edge detection');


%%
%STEP 5b: process detected edges 

%create a structure element
se1= strel('sphere',1); %Structure Element

%connect broken edges using MORPHOLOGICAL DILATION
result.edge=imdilate(result.edge,se1);
subplot(3,1,2);
imshow(result.edge);title('dilate to connect broken edges');

%fill up the inside of the edge aka the holes 
result.edge_filled = imfill(result.edge,4,'holes');
subplot(3,1,3);
imshow(result.edge_filled);title('fill up the inside of edges');


%%
%STEP 6: binarize the smoothened image (refer Step 4 )with a global threshold (NOT the edge image)

%compute a global threshold
result.threshold = graythresh(result.img_c);

%binarize
result.threshold_nuclei = imbinarize(result.img_c,result.threshold);

figure;
subplot(2,1,1);
imshow(result.threshold_nuclei);title('nuclei detected by thresholding');

%%
%STEP 7: find the nuclei that are NOT detected by thresholding

%create a BLACK balance image with the same size of the image 
result.balance_nuclei=false(size(result.edge_filled));

%compute the balance: 
%            balance nuclei = edge detected nuclei - thresholding detected nuclei
result.balance_nuclei=result.edge_filled-result.threshold_nuclei;
for row=1:v
    for col=1:h
        if result.balance_nuclei(row,col)<0 
            result.balance_nuclei(row,col)=0;
        end
    end
end

%convert balance from double -> logical for morpho processing later
%this is fine because edge_filled and threshold_nuclei are both logical
result.balance_nuclei=logical(result.balance_nuclei);

subplot(2,1,2);
imshow(result.balance_nuclei); title('nuclei not detected by thresholding');

%%
%STEP 8: processing the balance nuclei

figure; 
sgtitle('Morpho Processing nuclei undetected by thresholding');

% morphological processing
result.balance_nuclei=imerode(result.balance_nuclei,se1);
subplot(2,2,1);imshow(result.balance_nuclei); title('1.erode');

result.balance_nuclei=imopen(result.balance_nuclei,se1);
subplot(2,2,2);imshow(result.balance_nuclei); title('2.opening');

result.balance_nuclei=imopen(result.balance_nuclei,se1);
subplot(2,2,3);imshow(result.balance_nuclei); title('3.opening');

result.balance_nuclei=imfill(result.balance_nuclei,'holes');
subplot(2,2,4);imshow(result.balance_nuclei); title('4.fill holes');

%color label the nuclei
label_balance_nuclei = bwlabel(result.balance_nuclei);
result.color_balance_nuclei = label2rgb(label_balance_nuclei,@winter,'k','shuffle');

%%
%STEP 9: processing the thresholded nuclei

figure; 
sgtitle('Morpho processing nuclei detected by thresholding');

%separate connected nuclei 
%       subtracting the nuclei edge from the thresholded nuclei -> shrink the nuclei
result.threshold_nuclei = result.threshold_nuclei - result.edge;

[v,h] = size(result.threshold_nuclei);

for row=1:v
    for col=1:h
        if result.threshold_nuclei(row,col)<0 %negative
            result.threshold_nuclei(row,col)=0;
        end
    end
end

result.threshold_nuclei=logical(result.threshold_nuclei);
subplot(2,2,1);imshow(result.threshold_nuclei);title('1.edge subtraction');

result.threshold_nuclei=imerode(result.threshold_nuclei,se1);
subplot(2,2,2);imshow(result.threshold_nuclei);title('2. erode');

result.threshold_nuclei=imdilate(result.threshold_nuclei,se1);
subplot(2,2,3);imshow(result.threshold_nuclei); title('3. dilate');

result.threshold_nuclei=imfill(result.threshold_nuclei,'holes');
subplot(2,2,4);imshow(result.threshold_nuclei); title('4. fill holes');

%color label the nuclei
label_threshold_nuclei = bwlabel(result.threshold_nuclei);
result.color_threshold_nuclei = label2rgb(label_threshold_nuclei,@autumn,'k','shuffle');


%%
%display each color label of nuclei
figure;
subplot(2,1,1);
imshow(result.color_balance_nuclei);title('nuclei not detected by threshold in color');
subplot(2,1,2);
imshow(result.color_threshold_nuclei);title('nuclei detected by threshold in color');

%%
%STEP 10: combine both types of nuclei together as 1 image in each binary and color form
%           in color format, thresholding detected nuclei will have
%           one color map and thresholding UNdetected nuclei will have
%           another color map,  DISTINGUISH 
%       : find detected nuclei regions (connected components)
%       : find total count of nuclei regions (NumObjects)
    

result.combine_binary=result.threshold_nuclei + result.balance_nuclei; %BINARY FORMAT

for row=1:v
    for col=1:h
        if result.combine_binary(row,col)>0
            result.combine_binary(row,col)=1;
        else
            result.combine_binary(row,col)=0;
        end
    end
end

result.combine_color = result.color_threshold_nuclei + result.color_balance_nuclei; %COLOR FORMAT

%find detected nuclei regions
result.nuclei_components = bwconncomp(result.combine_binary,8);
%find total count of the detected nuclei regions
result.sum_nuclei=result.nuclei_components.NumObjects;

%display and compare before and after
figure;
imshowpair(imageinput,result.combine_binary,'montage');title('ORIGINAL vs DETECTED NUCLEI in BINARY');
figure;
imshowpair(result.combine_binary,result.combine_color,'montage'); title('ALL DETECTED NUCLEI: BINARY vs COLOR ');


%%
%FIND THE SIZE & SHAPE OF THE DETECTED NUCLEI and SHOW DISTRIBUTION

%get region properties / descriptors
nucleidata = regionprops(result.nuclei_components,'Area','Solidity','Eccentricity'); 

%%
%FIND SIZE, which is the area OR the number of pixels a region involves

%extract area data
result.nucleiarea = [nucleidata.Area]; 

%display in histogram : nucleus count against area of nucleus 
%                           number of nuclei having area x
%                           the distribution of each area involved
figure;
histogram(result.nucleiarea);
title('Size distribution of nuclei');
xlabel('area: number of pixels involved for each nucleus');
ylabel('nucleus count');


%%
%FIND SHAPE : shape descriptors


%extract solidity data (higher, more round), [0 1]
result.nucleisolidity = [nucleidata.Solidity]; 
%extract eccentricity data (higher, less round), [0 1]
result.nucleieccentricity = [nucleidata.Eccentricity];

%display solidity and circularity distribution
figure; sgtitle('Shape Descriptors and Distribution');
subplot(1,2,1);histogram(result.nucleisolidity);
xlabel('solidity');
ylabel('nucleus count');
subplot(1,2,2);histogram(result.nucleieccentricity);
xlabel('eccentricity');
ylabel('nucleus count');


%extract boundary coordinates of each nucleus region and plot
result.nuclei_bound = bwboundaries(result.combine_binary,8,'noholes');
[r,c]=size(result.nuclei_bound);
result.nuclei_chaincode=cell(r,c); %create empty cell for recording chaincode

figure; 
for i=1:length(result.nuclei_bound)
    
    %access the ith cell of nuclei_bound cell array, 
    %ith cell boundary data represents the boundary data of the ith nucleus
    bound=result.nuclei_bound{i};
    plot(bound(:,2),bound(:,1),'k'); %plot it horizontally plot(col,row)
    hold on;
    
    %generate chaincode of the nucleus boundary
    result.nuclei_chaincode{i} = fchcode(bound,8);%diagonal and non-diagonal Freemain chaincode
    %extract only the generated chaincode, ignore other data
    result.nuclei_chaincode{i} = result.nuclei_chaincode{i}.fcc; 
    
end
title('Shape: Boundary Plotting');
set(gca,'YDir','reverse'); %flipping plot at y-axis 


%%
%FIND THE BRIGHTNESS OF THE DETECTED NUCLEI AND SHOW DISTRIBUTION

%idea: 
%    1. nucleus_brightness : calculate the brightness of EACH nucleus by calculating the AVERAGE brightness of its IMAGE PIXELS 
%       because 1 nucleus may involve >1 image pixels  
%       and each image pixel may have different brightness value
%       so a good way to define the brightness of ONE nucleus is to find the average brightness of all its image pixels 
%
%    2. sum_nuclei_brightness : sum up the brightness of each nucleus computed above 
%       sum_nuclei_brightness += nucleus_brightness ;
%
%    3. avgbrightness : calculate the average brightness of ALL nuclei
%       avgbrightness = sum_nuclei_brightness / nuclei_count ;
%    

%et the total nuclei count
nuclei_count = result.sum_nuclei;

%make a tempo copy of the original GREEN EXTRACTED RGB image
result.imgtempo=result.img;

%initialize sum of brightness of all detected nuclei = 0 
result.sum_nuclei_brightness =0.00; 

%create an empty vector to store the brightness of each detected nucleus
%size will be the total nuclei count
result.nucleus_brightness = zeros(1,nuclei_count); 

%compute the brightness of all detected nuclei 1 by 1, -loop-
for nucleus_id=1:nuclei_count
    
    nucleus=false(size(result.combine_binary));%create a logical 2D black image
    nucleus(result.nuclei_components.PixelIdxList{nucleus_id}) = true; %show a detected nucleus region in binary ( based on nucleus_id )
    
    %extract that particular nucleus from the original RGB image and show it in the tempo RGB copy
    %tempo RGB copy will only have that ONE, CURRENTLY detected nucleus part shown in green, 
    %                              the surrounding will be black, including
    %                              the parts with other non-current nuclei
    for row=1:v
        for col=1:h
            if nucleus(row,col)== 1
                result.imgtempo(row,col,2) = result.img(row,col,2); 
            else
                result.imgtempo(row,col,2)=0;
            end
        end
    end
    
    %figure;imshowpair(nucleus,result.imgtempo,'montage');title('nucleus in binary mapped to temporary RGB image');
    
    %convert the tempo copy from RGB to HSV format
    %HSV : (:,:,1)=hue  (:,:,2)=saturation  (:,:,3)=brightness
    imgtempo_hsv= rgb2hsv(result.imgtempo);
    
    %extract the brightness data from the tempo HSV copy
    %this nucleus_pixels_brightness will be a [row,col] matrix  
    %                 recording the brightness value of ALL image pixels of the ENTIRE image, including 0 brightness pixels
    result.nucleus_pixels_brightness = imgtempo_hsv(:,:,3); 
    
    %since we only want to compute the brightness of the CURRENT NUCLEUS, 
    %we extract only the nucleus part which are those image pixels with a green hue brightness
    %the rest aka 0 brightness image pixels will be ignored
    %this nucleus_pixels_brightness will now be a vector
    %                   recording the brightness values of ONLY the image pixels involved in that one CURRENT NUCLEUS
    result.nucleus_pixels_brightness = result.nucleus_pixels_brightness(result.nucleus_pixels_brightness~=0); 
    
    %find the brightness of the current nucleus, 
    %       by calculating the average brightness of all its image pixels
    result.nucleus_brightness(1,nucleus_id) = sum(result.nucleus_pixels_brightness,'all') / numel(result.nucleus_pixels_brightness); 
    
    %append the sum of brightness of all detected nuclei, 1 by 1 each loop 
    result.sum_nuclei_brightness = result.sum_nuclei_brightness + result.nucleus_brightness(1,nucleus_id) ;
    
end

%compute final average brightness of all detected nucleii
%scalar value
result.nuclei_avgbrightness = result.sum_nuclei_brightness / nuclei_count;

%display in histogram -> [nucleus_count, brightness of nucleus] : the number of nuclei having that brightness
%                       the distribution of each invovled nucleus brightness
figure;
histogram(result.nucleus_brightness);
title('BRIGHTNESS distribution');
xlabel('brightness of nucleus'); xlim([0 1]);
ylabel('nucleus count');


end
    


    
                