%% 
figure
subplot(211)
m1=(mean(double(GL1),3)); % GL1 is the image
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
