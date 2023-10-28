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
r_data = zeros(image_h + 4, image_w + 4);
for k = 1:1:256
    r_data(2, k+2) = i_data(2, k);
    r_data(259, k+2) = i_data(255, k);
    r_data(k+2, 2) = i_data(k, 2);
    r_data(k+2, 259) = i_data(k, 255);
    r_data(2, 2) = i_data(2, 2);
    r_data(2, 259) = i_data(2, 255);
    r_data(259, 2) = i_data(255, 2);
    r_data(259, 259) = i_data(255, 255);

    r_data(1, k+2) = i_data(3, k);
    r_data(260, k+2) = i_data(254, k);
    r_data(k+2, 1) = i_data(k, 3);
    r_data(k+2, 260) = i_data(k, 254);
    r_data(1, 1) = i_data(3, 3);
    r_data(1, 260) = i_data(3, 254);
    r_data(260, 1) = i_data(254, 3);
    r_data(260, 260) = i_data(254, 254);

    r_data(1, 2) = i_data(3, 2);
    r_data(2, 1) = i_data(2, 3);
    r_data(1, 259) = i_data(3, 254);
    r_data(2, 260) = i_data(2, 254);
    r_data(259, 1) = i_data(255, 3);
    r_data(260, 2) = i_data(254, 2);
    r_data(259, 260) = i_data(255, 254);
    r_data(260, 259) = i_data(254, 255);
end
% disp(r_data(1, 1:258));
% disp(r_data);

% フィルタをかける用のデータ
m_data = zeros(image_h + 4, image_w + 4);
m_data = r_data;
m_data(3:(image_h + 2), 3:(image_w + 2)) = i_data;
% disp(m_data);

o_data = zeros(256);

% フィルタサイズ
filter_size = 5;

for x = 3:1:image_h+2
    for y = 3:1:image_w+2
        values = zeros(filter_size^2, 1);
        id = 1;
        for m = -2:1:2
            for n = -2:1:2
                i = x + m;
                j = y + n;

                values(id) = m_data(i, j);
                id = id + 1;
            end
        end
        % 中央値を計算
        median_value = median(values);
        o_data(x-2, y-2) = median_value;
        % disp(values);
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


% 以下はo_dataについての処理
disp("o_dataに関する値")
% mean
sum = 0;
for x = 1:image_h
    for y = 1:image_w
        sum = sum + o_data(x, y);
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
        tmp = (o_data(x, y) - mean)^2;
        sum = sum + tmp;
    end
end
variance = sum / (image_h * image_w);   
disp("variance");
disp(variance);

%フーリエスペクトルを計算
fs = fft2(o_data);
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
imshow(o_data,[gray_level_min gray_level_max]);
% disp(o_data)
%title('タイトル');


