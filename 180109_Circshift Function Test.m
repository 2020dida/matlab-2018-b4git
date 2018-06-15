clear all
%%

A=zeros([1024 50 5000]);

A(:,:,25)=ones([1024 50 1]);
Frame_Spacing_Lateral_Pixel_PreNPoint=1;

        for p=1:size(A,2)
                A(:,size(A,2)-p+1,:)=circshift(A(:,size(A,2)-p+1,:),Frame_Spacing_Lateral_Pixel_PreNPoint,3);
        end
        %%
        
        imagesc(A(:,:,25));
