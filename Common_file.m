clear;
clc;
close all;

img_name = '2';
img_path = 'imgs\2.jpg';
img = imread(img_path);
frac_dim_level_min = 2.3;
frac_dim_level_max = 3.3;
eig_cor_matrix_level = 1000;
level = 12;


DO_coord = task_1(img, img_name, frac_dim_level_min,frac_dim_level_max,eig_cor_matrix_level,level);
task_2(DO_coord,img, img_name )