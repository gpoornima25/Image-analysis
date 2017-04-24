%% 
figure
subplot(211)
% m=(double(GL1(:,:,1)));
% m1=mean(double(m),3);
m1=(mean(double(GL1),3));
imagesc(m1)

%Pixel shift correction (shift by one pixel for alternate odd rows)
subplot(212)
m2=m1;

nCol=size(m2,2)
nRow=size(m2,1)
for j=1:2:nRow
    for i=1: nCol-1
             m2(j,i)=m2(j,i+1);
    end
end
imagesc(m2)

 %%
%  %% 
% 
% subplot(311)
% m=(double(GL1(100:150,100:150,:)));
% m1=mean(double(m),3);
% imagesc(m1)
% 
% %shift by one pixel for alternate even rows
% subplot(312)
% m2=m1
% 
% n=size(m2);
% nCol=n(2)
% nRow=n(1)
% for j=2:2:nRow
%     for i=1: nCol-1
%              m2(j,i)=m2(j,i+1);
%     end
% end
% imagesc(m2)
% 
% %shift by 2 pixels for alternate even rows
% subplot(313)
% m2=m1
% 
% n=size(m2);
% nCol=n(2)
% nRow=n(1)
% for j=2:2:nRow
%     for i=1: nCol-2
%              m2(j,i)=m2(j,i+2);
%     end
% end
% imagesc(m2)