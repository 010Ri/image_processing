image_h=256; %画像の高さ
image_w=256; %画像の幅
gray_level_max= 255; %配布データの輝度値 (最大値)
gray_level_min= 0; %配布データの輝度値 (最小値)
input_file_name = 'idata1.dat'; %画像のファイル名
%ファイルオープン
fid=fopen(input_file_name,'r');
%画像データの読み込み
i_data = fread(fid,[image_w image_h],'uchar');
%今後，i_data を画像データの行列 (i_data(1,1) ～ i_data(256,256))として利用できる．
%縦横が逆転しているので，転置を取る．
i_data=i_data';
% disp(i_data(2, 1:256));
% disp(i_data(2, 1:256));
% disp(i_data);

% 折り返し用の行列の作成
r_data = zeros(image_h + 2, image_w + 2);
for k = 1:1:256
    r_data(1, k+1) = i_data(2, k);
    r_data(258, k+1) = i_data(255, k);
    r_data(k+1, 1) = i_data(k, 2);
    r_data(k+1, 258) = i_data(k, 255);
    r_data(1, 1) = i_data(2, 2);
    r_data(1, 258) = i_data(2, 255);
    r_data(258, 1) = i_data(255, 2);
    r_data(258, 258) = i_data(255, 255);
end
% disp(r_data(1, 1:258));
% disp(r_data);

% フィルタをかける用のデータ
m_data = zeros(image_h + 2, image_w + 2);
m_data = r_data;
m_data(2:(image_h + 1), 2:(image_w + 1)) = i_data;
% disp(m_data);

% Sobelフィルタ水平軸方向
sobel_h = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
% Sobelフィルタ垂直軸方向
sobel_v = [-1, -2, -1; 0, 0, 0; 1, 2, 1];

o_data_h = zeros(image_h, image_w); 
o_data_v = zeros(image_h, image_w);

for x = 2:1:image_h+1
    for y = 2:1:image_w+1
        sum_h = 0;
        sum_v = 0;
        for m = -1:1:1
            for n = -1:1:1
                i = x + m;
                j = y + n;
                sum_h = sum_h + m_data(i, j) * sobel_h(m + 2, n + 2);
                sum_v = sum_v + m_data(i, j) * sobel_v(m + 2, n + 2);
            end
        end
        o_data_h(x - 1, y - 1) = sum_h;
        o_data_v(x - 1, y - 1) = sum_v;
    end
end

% 以下はi_dataについての処理
disp("i_dataに関する値")
% mean
sum = 0;
for x = 1:image_h
    for y = 1:image_w
        sum = sum + i_data(x, y);
    end
end
mean = sum / (image_h * image_w);   
disp("mean");
disp(mean);

% variance
sum = 0;
tmp = 0;
for x = 1:image_h
    for y = 1:image_w
        tmp = (i_data(x, y) - mean)^2;
        sum = sum + tmp;
    end
end
variance = sum / (image_h * image_w);   
disp("variance");
disp(variance);

% シグマの計算
sigma = sqrt(variance);
disp("sigma")
disp(sigma);
t_sigma = 3 * sigma;
disp("3sigma")
disp(t_sigma);

%フーリエスペクトルを計算
fs = fft2(i_data);
%エネルギースペクトルを計算
ps = fs.*conj(fs);
% auto-correlation function
ac = ifft2(ps);
ac_i = real(ac)/(image_h*image_w);
c_f_i = ac_i - (mean^2);
% auto-correlation coefficient
acc_i = c_f_i / variance;

% 最大値と最小値
min_h = min(o_data_h, [], "all");
max_h = max(o_data_h, [], "all");
min_v = min(o_data_v, [], "all");
max_v = max(o_data_h, [], "all");

% 輝度補正
A = min_h;
B = max_h;
% C = 0;
% D = 255;
C = mean - t_sigma;
D = mean + t_sigma;
for x = 1:256
    for y = 1:256
        o_data_h(x,y) = (D - C) / (B - A) * (o_data_h(x,y) - A) + C;
    end
end

A = min_v;
B = max_v;
% C = 0;
% D = 255;
C = mean - t_sigma;
D = mean + t_sigma;
for x = 1:256
    for y = 1:256
        o_data_v(x,y) = (D - C) / (B - A) * (o_data_v(x,y) - A) + C;
    end
end

% 以下はo_dataについての処理
disp("o_dataに関する値")
% mean
sum = 0;
for x = 1:image_h
    for y = 1:image_w
        sum = sum + o_data_h(x, y);
    end
end
mean = sum / (image_h * image_w);   
disp("mean");
disp(mean);

% variance
sum = 0;
tmp = 0;
for x = 1:image_h
    for y = 1:image_w
        tmp = (o_data_h(x, y) - mean)^2;
        sum = sum + tmp;
    end
end
variance = sum / (image_h * image_w);  
disp("variance");
disp(variance);

%フーリエスペクトルを計算
fs = fft2(o_data_h);
%エネルギースペクトルを計算
ps = fs.*conj(fs);
% auto-correlation function
ac = ifft2(ps);
ac_o = real(ac)/(image_h*image_w);
c_f_o = ac_o - (mean^2);
% auto-correlation coefficient
acc_o = c_f_o / variance;

dif = acc_o - acc_i;
% disp(dif); 


%画像データの表示 -- 原画像の表示
imshow(o_data_h,[gray_level_min gray_level_max]);
% imshow(o_data_v,[gray_level_min gray_level_max]);
% disp(o_data)
%title('タイトル');


