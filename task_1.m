clear;
clc;
close all;

img_path = 'imgs\2.jpg';
data_save_name = '2_img_coeffs'
img = imread(img_path);
frac_dim_level = 2.5;
eig_cor_matrix = 2500;
level = 10;
bin_level = 0.3;


img = rgb2gray(img);
img = im2double(img);

img = medfilt2(img);

v_k_img_vertical = img;
v_k_img_gorizontal = img;
v_k_img = img;

[M,N] = size(img);

for i = 1:M
    v_k_img_gorizontal(i,:) = cwt(img(i,:),8, 'haar');
end

for i = 1:N
    v_k_img_vertical(:,i) = cwt(img(:,i),8, 'haar');
end

for i = 1:M
    for j = 1:N
        if(v_k_img_gorizontal(i, j) ~= v_k_img_vertical(i, j))
            if(abs(v_k_img_gorizontal(i, j)) > abs(v_k_img_vertical(i, j)))
                v_k_img(i, j) = v_k_img_vertical(i, j);
            else
                v_k_img(i, j) = v_k_img_gorizontal(i, j);
            end
        else
            v_k_img(i, j) = v_k_img_vertical(i, j);
        end
    end
end

mass_img_none_all = img;

for i = 1:M
    for j = 1:N
        if(abs(v_k_img(i, j))*255 < level)
            mass_img_none_all(i, j) = mass_img_none_all(i, j)*0;
        end
    end
end

BW = im2bw(mass_img_none_all, bin_level);
imshow(BW)
BW = imfill(BW, 8, 'holes');

for k = 1:50
    for i = 1:M-1
        for j=1:N-1
            if((BW(i+1, j+1) == 1) & (BW(i, j) == 1))
                BW(i, j+1) = 1;
                BW(i+1, j) = 1;
            end
            if((BW(i, j+1) == 1) & (BW(i+1, j) == 1))
                 BW(i+1, j+1) = 1;
                 BW(i, j) = 1;
            end
        end
    end
end

BW2 = BW;
v = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];

for i = 1:M
    for j = 1:N
        if(BW2(i, j) == 1)
            k = 1;
            while (((i+k) <= M) & (BW2(i + k, j) == 1))
                k = k + 1;
            end
            d = 1;
            while(((j+d) <= N) & (BW2(i, j+d) == 1))
                d = d + 1;
            end
            if((d ~= 1) & (k ~= 1))
                BW2(i:i + k, j:j + d) = 0;
                v = [v;i, i+k, j, j+d, 0, 0, 0, 0, 0, 0;];
            end
        end
    end
end


[A, Y] = size(v);
output = img;
for i = 2:A   
    output = insertShape(output, 'rectangle',[v(i,3), v(i,1), v(i,4)-v(i,3), v(i,2)-v(i,1)] ,'LineWidth',2);
end
imshow(output)

for i = 2:A
    mass = img(v(i,1):v(i,2)-1, v(i,3):v(i,4)-1) * 255;
    [M,N] = size(mass);
    if (M > N) 
        a = N;
    else
        a = M;
    end
    for Q = 1:a
        tx = floor(M/Q);
        ty = floor(N/Q);
        sum = 0;
        for X = 1:tx
            for Y = 1:ty
                m = min(min(mass((Q*(X-1)+1):(Q*X), (Q*(Y-1)+1):(Q*Y))));
                sum = sum + ceil((m)/(Q));
            end
        end
        if (sum == 0)
            sum = 1;
        end
        sumfr(Q) = sum;
    end
    A11 = a; 
    A12 = 0;
    A21 = 0;
    A22 = 0;
    d1 = 0;
    d2 = 0;
    for iv = 1:a
        A12 = A12 + log(iv);
        A22 = A22 + log(iv)*log(iv);
        d1 = d1 + log(sumfr(iv));
        d2 = d2 + log(sumfr(iv))*log(iv);
    end
    A21 = A12;
    D = A11*A22 - A21*A12;
    D1 = d1*A22 - d2*A12;
    D2 = A11*d2 - A21*d1;
    b = D1/D;
    m = D2/(d*1.4);
    v(i, 6) = abs(m);
end

for i = 2:A
    cam = img(v(i,1):v(i,2)-1, v(i,3):v(i,4)-1) * 255;
    [M,N] = size(cam);
    if(M>N) N=M; end
    if(M<N) M=N; end
    sigma2_x = var(cam(:));
    mean_x = mean(cam(:));
    cam_r = circshift(cam,[1 0]);
    cam_c = circshift(cam,[0 1]);
    rho_mat = corrcoef([cam(:);cam(:)],[cam_c(:); cam_r(:)]);
    rho = rho_mat(1, 2);
    [rr, cc] = meshgrid([-M:M-1], [-N:N-1]);
    r_x = (sigma2_x * rho.^sqrt(rr.^2+cc.^2));
    mz = max(eig(real(r_x)));
    v(i,5) = mz;
end
e = [];
for i = 2:A
    if((v(i,5) >= eig_cor_matrix) & (v(i,6) >= frac_dim_level))
        r = 1;
    else
        r = 0;
    end
    e = [e;i,r;];
end
DO = [];
DO_coord = [];
for i = 1:A
    if (e(i, 2) == 1)
        DO = [DO; v(i,:)];
        DO_coord = [DO_coord;  v(i,3), v(i,1), v(i,4)-v(i,3), v(i,2)-v(i,1);];
    end
end
    

output = insertShape(img, 'rectangle',DO_coord ,'LineWidth',2);
imshow(output)



