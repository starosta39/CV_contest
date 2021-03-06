function [DO_coord] = task_1(img, img_name, frac_dim_level_min,frac_dim_level_max,eig_cor_matrix_level,level )
img_col = img;
img = rgb2gray(img);
img = im2double(img);

img = medfilt2(img);
mean_bright = mean2(img);
img_std = std2(img);
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

BW = imbinarize(mass_img_none_all,'adaptive');
BW = imfill(BW, 8, 'holes');
for k = 1:50
    for i = 1:M-1
        for j=1:N-1
            if((BW(i+1, j+1) == 1) && (BW(i, j) == 1))
                BW(i, j+1) = 1;
                BW(i+1, j) = 1;
            end
            if((BW(i, j+1) == 1) && (BW(i+1, j) == 1))
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
            while (((i+k) <= M) && (BW2(i + k, j) == 1))
                k = k + 1;
            end
            d = 1;
            while(((j+d) <= N) && (BW2(i, j+d) == 1))
                d = d + 1;
            end
            if((d ~= 1) && (k ~= 1))
                BW2(i:i + k, j:j + d) = 0;
                v = [v;i, i+k-1, j, j+d-1, 0, 0, 0, 0, 0, 0;];
            end
        end
    end
end

[A, L] = size(v);

for i = 2:A
    mass = img(v(i,1):v(i,2), v(i,3):v(i,4)) * 255;
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
    m = D2/(D*1.4);
    v(i, 6) = abs(m);
end

for i = 2:A
    cam = img(v(i,1):v(i,2), v(i,3):v(i,4)) * 255;
    [M,N] = size(cam);
    if(M>N) N=M; end
    if(M<N) M=N; end
    sigma2_x = var(cam(:));
    mean_x = mean(cam(:));
    cam_r = circshift(cam,[1 0]);
    cam_c = circshift(cam,[0 1]);
    rho_mat = corrcoef([cam(:);cam(:)],[cam_c(:); cam_r(:)]);
    rho = rho_mat(1, 2);
    if isnan(rho)
        rho = 0;
    end
    
    [rr, cc] = meshgrid([-M:M-1], [-N:N-1]);
    r_x = (sigma2_x * rho.^sqrt(rr.^2+cc.^2));
    mz = max(eig(real(r_x)));
    v(i,5) = mz;
    v(i,7) = mean_x;
    v(i,8) = sigma2_x;
end



e = [0,0];
for i = 2:A
    if(v(i,5) >= eig_cor_matrix_level && (v(i,6) >= frac_dim_level_min && v(i,6) <= frac_dim_level_max))
        r = 1;
    else
        r = 0;
    end
    e = [e;i,r;];
end
DO = [0,0,0,0,0,0,0,0,0,0];
DO_coord = [0,0,0,0];
for i = 1:A
    if (e(i, 2) == 1)
        DO = [DO; v(i,:)];
        DO_coord = [DO_coord;  v(i,3), v(i,1), v(i,4)-v(i,3), v(i,2)-v(i,1);];
    end
end

output = img_col;
for i = 2:A   
    output = insertShape(output, 'rectangle',[v(i,3), v(i,1), v(i,4)-v(i,3), v(i,2)-v(i,1)] ,'LineWidth',2);
    %output = insertText(output,[v(i,3),v(i,1)], num2str(i),'FontSize', 10, 'TextColor','green','BoxOpacity',0);
end
figure,imshow(output);
all_DO = output;

[n_real_DO, n_marks] = size(DO_coord);
output = img_col;
for i = 2:n_real_DO  
    output = insertShape(output, 'rectangle',DO_coord, 'LineWidth',2);
    %output = insertText(output,[DO_coord(i,1),DO_coord(i,2)], num2str(i),'FontSize', 10, 'TextColor','green','BoxOpacity',0);
end
figure,imshow(output);
true_DO = output;

answ = input('Do you want enter data in the table ?  1-YES/0-NO ');
if (answ == 1)
    if (exist('True_DO_coefs.csv', 'file') == 2)&&(exist('All_DO_coefs.csv', 'file') == 2)
        True_DO_coefs = readtable('True_DO_coefs.csv');
        All_DO_coefs = readtable('All_DO_coefs.csv');        
    else
        Number_on_img = []; 
        Img_name = []; 
        X = [];
        Y = [];
        Width = [];
        High  = [];
        Window_mean = [];
        Window_std = [];
        Img_mean = [];
        Img_std = [];
        Wc_coef_level = [];
        Fractal_dim_level_min = [];
        Fractal_dim_level_max = [];
        Eigenvalue_of_matrix_level = [];
        Fractal_dim = [];
        Eigenvalue_of_matrix = [];
        Answer = [];
        True_DO_coefs = table(Number_on_img,Img_name,X,Y,Width, High, Window_mean, Window_std, Img_mean, Img_std, Wc_coef_level,Fractal_dim_level_min,Fractal_dim_level_max,Eigenvalue_of_matrix_level,Fractal_dim,Eigenvalue_of_matrix,Answer);              
        All_DO_coefs = True_DO_coefs;                        
    end
    for i = 2:A 
        cell_DO = {i-1, img_name, v(i,3), v(i,1), v(i,4)-v(i,3), v(i,2)-v(i,1),v(i,7), v(i,8), mean_bright, img_std, level, frac_dim_level_min,frac_dim_level_max, eig_cor_matrix_level, v(i,6),v(i,5), e(i, 2)};
        All_DO_coefs = [All_DO_coefs; cell_DO];
        if (e(i, 2) == 1)
            True_DO_coefs = [True_DO_coefs; cell_DO];
        end
    end
    close('all')
    writetable(True_DO_coefs,'True_DO_coefs.csv');
    close('all')
    writetable(All_DO_coefs, 'All_DO_coefs.csv');
    close('all')
end
answ = input('Do you want save imgs ?  1-YES/0-NO ');
if (answ == 1)
    if (exist('All_DO', 'dir') == 7 && exist('True_DO', 'dir') == 7)
        imwrite(all_DO ,strcat('All_DO\',img_name,'all_detected_DO.png'))
        imwrite(true_DO ,strcat('True_DO\',img_name,'true_detected_DO.png'))
    else 
        mkdir All_DO;
        mkdir True_DO;
        imwrite(all_DO ,strcat('All_DO\',img_name,'all_detected_DO.png'))
        imwrite(true_DO ,strcat('True_DO\',img_name,'true_detected_DO.png'))
    end
end
end
