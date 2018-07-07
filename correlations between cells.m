

%% correct the shift 
for j = 1:H2d/2
    Data2Dshift(j*2,2:W2d,:) = Data2D(j*2,1:W2d+1-2,:);
end
Ca2dImageShift = mean(Data2Dshift(:,:,:),3);
figure; imagesc(Ca2dImageShift);title('shift');

%% detect signle cell

sm = min(min(Ca2dImage));
lg = max(max(Ca2dImage));
r = lg-sm;
mask = Ca2dImage > r/10;
% Ca2d = uint8(round((Ca2dImage-sm)/r*255));
Ca2d = Ca2dImage;
background6 = imopen(Ca2d,strel('disk',2));
% background6 = imopen(Ca2d,strel('disk',15));
Ca2da6 = imadjust(Ca2d-background6);
background5 = imopen(Ca2da6,strel('disk',13));
Ca2da5 = imadjust(Ca2da6-background5);
background4 = imopen(Ca2da5,strel('disk',11));
Ca2da4 = imadjust(Ca2da5-background4);
background3 = imopen(Ca2da4,strel('disk',9));
Ca2da3 = imadjust(Ca2da4-background3);

figure
subplot(231); imagesc(Ca2d);title('original');
subplot(232); imagesc(Ca2da6); title('disk 15');
subplot(233); imagesc(Ca2da5); title('disk 13');
subplot(234); imagesc(Ca2da4); title('disk 11');
subplot(235); imagesc(Ca2da3); title('disk 9');

figure
subplot(231); imagesc(Ca2d);title('original');
subplot(232); imagesc(background6); title('b disk 15');
subplot(233); imagesc(background5); title('b disk 13');
subplot(234); imagesc(background4); title('b disk 11');
subplot(235); imagesc(background3); title('b disk 9');


Ca2db = Ca2da6;
% Ca2db = redfigure;
% H =  fspecial('disk',2);
H = [0 1 0; 1 1 1; 0 1 0]/5;
Ca2dM = uint8(zeros(size(Ca2d)));
Cells = [];
% detect cells
for i = 1:256*0.8
    level = (255-i)/255;
    Ca2dc = im2bw(Ca2db,level);
    Ca2dd = bwareaopen(Ca2dc, 30); %remove <30 pixel detections
    Ca2dM = Ca2dM+uint8(Ca2dd);   % 
    % Ca23d = imfill(Ca23d,'holes');
    %     subplot(6,8,i);
    %     imagesc(Ca23d);
    cc = bwconncomp(Ca2dd, 8);
    
    num = cc.NumObjects;
    
    ccsingle.Connectivity = cc.Connectivity;
    ccsingle.ImageSize = cc.ImageSize;
    ccsingle.NumObjects = 1;
    
    if num >0
        % dd = regionprops(cc, 'basic');
        for j = 1:num
            % dd(j).PixelIdxList = cc.PixelIdxList(j);
            pixels = cc.PixelIdxList{j};
            cellarea = length(pixels);
            
            if cellarea > 49
                w = zeros(H2d,W2d);
                w(pixels) = 1;
                w = imfill(w);
                pixelsNew = find(w == 1);
                ccsingle.PixelIdxList{1} = pixelsNew;
                dd = regionprops(ccsingle, 'basic');
                dd.PixelIdxList = ccsingle.PixelIdxList;
                Cells = [Cells dd];
            else
                Ca2dd(pixels) = 0;
            end
        end
        
    end
    
    Ca23d1 = imfilter(double(Ca2dd),H);
    % Ca23d1 = imfilter(double(Ca2dd),H);
    % Ca23d2 = imfilter(double(Ca23d),H1);
    Ca23d3 = uint8(Ca23d1==0);
    % Ca23b = uint8(double(Ca23b.*Ca23d3).*(1-Ca23d1));
    Ca2db = Ca2db.*double(Ca23d3);
end
length(Cells)

% % check if the cells detected are in cell like shape 
%  
% ind  = 1
% while ind < length(Cells);
%     c1 = Cells(ind).Centroid;
%     pixels1 = Cells(ind).PixelIdxList{:};
%     w = zeros(H2d,W2d


% check if there will be overlapped cells, based on distance
ind = 1;
while ind <length(Cells)
    c1 = Cells(ind).Centroid;
    com = ind+1;
    while com <=length(Cells)
        c2 = Cells(com).Centroid;
        d = sqrt((c1(1)-c2(1))^2+(c1(2)-c2(2))^2); %distance(c1,c2);
        if d < 7.5 
            disp([num2str(d) '--' num2str(ind) '-' num2str(com)])
            if com == length(Cells)
                Cells = Cells(1:com-1);
            else
                Cells = [Cells(1:com-1) Cells(com+1:length(Cells))];
            end
            ind =1;
            c1 = Cells(ind).Centroid;
            com = ind;
        end
        com = com +1;
    end
    ind = ind +1;
end= ind +1;
% end

%% pull out single cell trace

% if W2d < 256
%     Data2D = imresize(Data2D,2);
% end

clear CellTrace CellPosition;

Data2 = reshape(Data2D,[H2d*W2d frames]);
% Data2 = reshape(a1,[256*256 frames]);
for i = 1:length(Cells)
    a = Cells(i).PixelIdxList{:};
    % a = fix((a-1)/4)+1;
    b = sum(Data2(a,:))/Cells(i).Area;
    Cells(i).Trace = b;
    CellPosition(:,i) = Cells(i).Centroid;
    CellTrace(:,i) = b;
end


l = length(Ca2Dfile);
Sname = [Ca2Dfile(1:l-3) '_Cells.mat'];
save(Sname,'Cells');

h= figure; imagesc(Ca2d); hold on;
plot(CellPosition(1,:),CellPosition(2,:),'ro','linewidth',1);
saveas(h,[Ca2Dfile(1:l-3) 'fig'],'fig');

%savefig(h,[Ca2Dfile(1:l-3) 'fig']);

CellTrace_backup = CellTrace;

% CellTrace = CellTrace_backup; %(228:298,:);
[R P] = corrcoef(CellTrace);

for i = 1:length(Cells)
    R(i,i) = 0;
    P(i,i) = 1;
end

h=figure;
imagesc(R); colorbar; title('R');
% saveas(h,[Ca2Dfile(1:l-4) '_R.fig'],'fig');
CellTrace = CellTrace_backup;

% for c = 1:1 %length(LineCellPosition)
%     COIPositionY = yst + LineCellPosition(c);
%     clear 'CellDistance';
%     for i = 1:length(CellArea)
%         CellDistance(i) = distance(CellPosition(i,1),CellPosition(i,2),COIPositionX,COIPositionY);
%     end
%     [d CellInd] = min(CellDistance);
% end
%
% %%
% for i = 1:543
%     xst = max([1 fix(CellPosition(i,1)-1)]);
%     xed = min([128 fix(CellPosition(i,1)+1)]);
%     yst = max([1 fix(CellPosition(i,2)-1)]);
%     yed = min([128 fix(CellPosition(i,2)+1)]);
%     trace = squeeze(mean(mean(a1(xst:xed,yst:yed,:))));
%     CellTrace(:,i) = trace;
% end
%
% for i = 1:542
%     Trace1 = CellTrace(:,i);
%     for j = i+1:543
%         Trace2 = CellTrace(:,j);
%         [c,lags] = xcorr(Trace1,Trace2);
%         [M I] = max(c);
%         rmax(i,j) = M; dmax(i,j) = lags(I);
%         [M I] = min(c);
%         rmin(i,j) = M; dmin(i,j) = lags(I);
%     end
% end
%



Fs = 30;
[n_T,n_channel] = size(CellTrace);
% Convert it into delta(f)/f
T_baseline = 5;
f_baseline = mean(CellTrace(1:T_baseline*Fs,:),1);
S = CellTrace ./ repmat(f_baseline,[n_T 1]);
% plot(S); % delta(f)/f
 

% Filter out the high frequency noise
band_pass = 2;
Filter_order = 512;
[b,a] = fir1(Filter_order,band_pass/(Fs/2));
S = filtfilt(b,a,S);

wo = 1.25/(Fs/2);  bw = wo;
[b,a] = iirnotch(wo,bw);
S = filtfilt(b,a,S);

[l n] = size(CellTrace); clear PeakValue;

for i =  1:n
    a = cwt(S(:,i),100,'haar');
    PeakValue(i) = a(580)-a(652);
end

[q w] = sort(PeakValue);

figure
for i = 1:20
    subplot(5,4,i)
    plot(CellTrace(:,w(i)),'r');
    title(num2str(w(i)));
end


figure
for i = 1:20
    subplot(5,4,i)
    plot(CellTrace(:,w(n-20+i)),'k');
    title(num2str(w(n-20+i)));
end

figure
imagesc(Ca2d);
hold on
plot(CellPosition(1,w(1:20)),CellPosition(2,w(1:20)),'r*');
plot(CellPosition(1,w(n-19:n)),CellPosition(2,w(n-19:n)),'wo');
hold off;

figure
imagesc(ColorImage);
hold on
plot(CellPosition(1,w(1:20)),CellPosition(2,w(1:20)),'r*');
plot(CellPosition(1,w(n-19:n)),CellPosition(2,w(n-19:n)),'wo');
hold off;