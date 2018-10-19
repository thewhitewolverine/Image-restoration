clc;
close all;
clear all;
img=im2double(imread('Blurry.jpg'));
blur=im2double(rgb2gray(imread('Blurr1_1_result_Cho_5.png.psf.png')));
im_r=img(:,:,1);
im_g=img(:,:,2);
im_b=img(:,:,3);
imgray=rgb2gray(img);
[M,N]=size(imgray);
u=0:1:M-1;
v=0:1:N-1;
W_N=exp(-1j*2*3.1416/N);
W_M=exp(-1j*2*3.1416/M);
Wm=u.^W_M;
Wn=v.^W_N;
ft=imgray;
bft=[];
[Mb,Nb]=size(blur);


for i=1:M
    for j=1:N
        bft(i,j)=0;
        if(i<=Mb)
            if(j<=Nb)
                bft(i,j)=blur(i,j);
            end
        end
    end
end
bftf=bft;
for i=1:M
    for j=1:N
        ft(i,j)=0;
        bft(i,j)=0;
        sum=0.000;
        for x=1:M
            for y=1:N
                sum=sum+double(imgray(x,y))*(W_M^((i-1)*(x-1)))*(W_N^((j-1)*(y-1)));
                
            end
        end
        ft(i,j)=sum;
    end
end


