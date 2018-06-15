clear all
%%

N=4;

S=(N)^0.5;

Image_Input=imread('C:\TuanShu\TestImage.jpg');

Image_LAB=rgb2lab(Image_Input);

X_Spacing=size(Image_LAB,1)/S;
Y_Spacing=size(Image_LAB,2)/S;
C=zeros(5,N);   %l a b x y

for p=1:N
   C(1:3,p)=Image_LAB
end


