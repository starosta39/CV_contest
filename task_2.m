function [] = task_2(DO_coord, img, img_name)
img_col = img;
img = rgb2gray(img);
img = im2double(img);
comm_mask = zeros(size(img));
[n_real_DO, n_marks] = size(DO_coord);
for i  = 2:n_real_DO
    I = img(DO_coord(i,2):DO_coord(i,2) + DO_coord(i,4),DO_coord(i,1):DO_coord(i,1) + DO_coord(i,3));
    [M,N] = size(I);
    %I = imgaussfilt(I,1);
    %I = medfilt2(I);
    hy = fspecial('sobel');
    hx = hy';
    Iy = imfilter(double(I), hy, 'replicate');
    Ix = imfilter(double(I), hx, 'replicate');
    gradmag = sqrt(Ix.^2 + Iy.^2);
    BW_DO = edge(gradmag, 'canny');   
    se90 = strel('line',3,90);
    se0 = strel('line',3,0);
    BWsdil = imdilate(BW_DO,[se90 se0]);
    BWdfill = imfill(BWsdil, 'holes');
    new_field = zeros(max([M,N])*2);
    start_i = int32(size(new_field)/4);
    new_field(start_i(1):start_i(1) + M -1 ,start_i(1):start_i(1) + N - 1) = BWdfill;
    BWnobord = imclearborder(new_field,4);
    seD = strel('diamond',1);
    BWfinal = imerode(BWnobord,seD);
    BWfinal = imerode(BWfinal,seD);
    mask = new_field(start_i(1):start_i(1) + M -1, start_i(1):start_i(1) + N -1  );
    comm_mask(DO_coord(i,2):DO_coord(i,2) + DO_coord(i,4),DO_coord(i,1):DO_coord(i,1) + DO_coord(i,3)) = mask;
    
end
res = labeloverlay(img_col,comm_mask,'Colormap',[1,1,0],'Transparency',0.6);
figure, imshow(comm_mask);
figure, imshow(res);

answ = input('Do you want save imgs ?  1-YES/0-NO ');
if (answ == 1)
    if (exist('DO_zone', 'dir') == 7 )
        imwrite(res ,strcat('DO_zone\',img_name,'true_detected_DO_zone.png'))
    else 
        mkdir DO_zone;
        imwrite(res ,strcat('DO_zone\',img_name,'true_detected_DO_zone.png'))
    end
end

