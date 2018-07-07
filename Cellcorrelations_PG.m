% Hongtao Ma,  detected the cells location and calculate the Corrcoef among
% Cells.

[Ca2Dfile, pathname] = uigetfile('*.blk', 'Please select the 2D scan file');
cd(pathname);
[Data2D fileinfo] =  fileopen(Ca2Dfile);
[H2d W2d frames] = size(Data2D);

Ca2dImage = mean(Data2D(:,:,:),3);
figure;imagesc(Ca2dImage);title('original');

%% correct the shift 
for j = 1:H2d/2
    Data2Dshift(j*2,2:W2d,:) = Data2D(j*2,1:W2d+1-2,:);
end
Ca2dImageShift = mean(Data2Dshift(:,:,:),3);
figure; imagesc(Ca2dImageShift);title('shift');
%%

% Ca2dImage = mean(Data2D(:,:,500:2000),3);

% [CaLineScanFile pathname] = uigetfile('*.mat', 'Please select the line scan file');
% load(CaLineScanFile);
%
% [CaLineScanImage pathname] = uigetfile('*.tif', 'Please select the line scan image');
%
% CaLSImage = importdata([pathname CaLineScanImage]);
% % CaLSImage = imrotate(CaLSImage,-90);
% CaLSImage = imresize(CaLSImage,4/5);
% CaLSImage = double(CaLSImage);
%
% %% find the overlab between two file
%
% % find the poential best correlation part
% c = xcorr2(double(Ca2dImage), double(CaLSImage));
% [i,j] = find(c > max(max(c))*0.75);
% s1 = min(i); s2 = min(j);
% e1 = max(i); e2 = max(j);
%
% % find the best overlaping coordinate
% WLS = length(CaLSImage);
% [l w] = size(c);
% co = zeros(size(c));
% for i = s1:e1
%     Xst2d = max([1 i-WLS+1]);
%     Xed2d = min([i W2d]);
%     XstLS = max([1 WLS-i+1]);
%     XedLS = min([l-i+1 WLS]);
%     for j = s2:e2
%         Yst2d = max([1 j-WLS+1]);
%         Yed2d = min([j W2d]);
%         YstLS = max([1 WLS-j+1]);
%         YedLS = min([l-j+1 WLS]);
%         S2d = double(Ca2dImage(Xst2d:Xed2d,Yst2d:Yed2d));
%         SLS = CaLSImage(XstLS:XedLS,YstLS:YedLS);
%         q = corrcoef(S2d,SLS);
%         co(i,j) = q(1,2);
%     end
% end
%
% [mx2 i] = max3(co);
% xc = i(2); yc = i(1);
%
%     Xst2d = max([1 xc-WLS+1]);
%     Xed2d = min([xc W2d]);
%     Yst2d = max([1 yc-WLS+1]);
%     Yed2d = min([yc W2d]);
% CutCaImage = Ca2dImage(Yst2d:Yed2d,Xst2d:Xed2d);


% shiftx = l1-(xc-(xLS-1)); shifty = yc-(yLS-1);
% xst = fix((l1-xLS)/2)-shiftx; yst = fix((w1-yLS)/2) + shifty;

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



% Ca2db = imadjust(background4);
% Ca2db = Ca2da3;
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
%             ShowSingleCell(ind,Cells,Ca2d);
%             ShowSingleCell(com,Cells,Ca2d);
%             pixels1 = Cells(ind).PixelIdxList{:};
%             pixels2 = Cells(com).PixelIdxList{:};
%             w = zeros(H2d,W2d);
%             w(pixels1) = 1; w(pixels2) = 1; w = imfill(w);
%             pixelsNew = find(w == 1);
%             ccsingle.PixelIdxList{1} = pixelsNew;
%             dd = regionprops(ccsingle, 'basic');
%             dd.PixelIdxList = ccsingle.PixelIdxList;
%             Cells(ind) = dd;
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
end

% while Cells(1).Area <10
%     l = length(Cells);
%     Cells = Cells(2:l);
% end
%
%
% ind = 2;
% while ind < length(Cells)
%     if Cells(ind).Area <10
%         % disp(num2str(ind));
%         l = length(Cells);
%         Cells= [Cells(1:ind-1) Cells(ind+1:l)];
%         ind = ind -1;
%     end
%     ind = ind +1;
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
% for lev = 1:1
%     le = max(max(R))*0.1*(10-lev);
%     % le = max(max(R))*0.1*6.5;
%     R1 = R>le;
%     q = find(R1>0);
%     w = P(q);
%     if max(w) < 0.01        
%         [I,J] = find(R1 >0);        
%         group  = 1;
%         ind = 1;
%         cellpair = 1;
%         clear CellGroup;
%         while cellpair <= length(I)
%             clear Gcell;
%             Gcell(1) = I(cellpair);
%             Gcell(2) = J(cellpair);
%             I(cellpair) = 0;
%             J(cellpair) = 0;
%             Gcell = sort(Gcell);
%             ind = 1;
%             if Gcell(1) > 0
%                 while ind <= length(Gcell)
%                     l = length(Gcell);
%                     k = Gcell(ind);
%                     g1 = find(I == k);
%                     if ~isempty(g1)
%                         for i = 1:length(g1)
%                             Gcell(i+l) = J(g1(i));
%                             I(g1(i)) = 0;
%                             J(g1(i)) = 0;
%                         end
%                     end
%                     l = length(Gcell);
%                     
%                     g2 = find(J == k);
%                     if ~isempty(g2)
%                         for i = 1:length(g2)
%                             Gcell(i+l) = I(g2(i));
%                             I(g2(i)) = 0;
%                             J(g2(i)) = 0;
%                         end
%                     end
%                     
%                     if isempty([g1' g2'])
%                         ind = ind + 1;
%                         %                 else
%                         %                     Gcell = sort(Gcell);
%                         %                     ind = 1;
%                     end
%                 end
%                 Gcell = unique(Gcell);
%                 CellGroup{group} = Gcell;
%                 group = group + 1;
%                 ind = 1;
%             end
%             cellpair = cellpair + 1;
%         end
%         
%         h = figure;
%         imagesc(Ca2d); colormap('gray'); hold on;
%         
%         ColorIndex = 'rgbcmykwrgbcmykw';
%         cellshape ='+osd^v<>ph';
%         for i = 1:length(CellGroup)
%             cellind = CellGroup{i};
%             q1 = [ColorIndex(mod(i,8)+1) cellshape(mod(i,10)+1)];
%             plot(CellPosition(1,cellind),CellPosition(2,cellind),q1,'linewidth',2);
%         end
%         title([num2str(0.1*(9-lev)) ' max = ' num2str(le)]);
%         l = length(Ca2Dfile);
%         saveas(h,[Ca2Dfile(1:l-4) '_' num2str(le) '.fig'],'fig');
%         h = figure;
%         [Frames n] = size(CellTrace);
% 
%         x = 1:Frames;
%         % x = x*0.256;
%         for i = 1:length(CellGroup)
%             t = zeros(size(CellTrace_backup(:,1)));            
%             subplot(ceil(length(CellGroup)/2),2,i); hold on;
%             q1 = [ColorIndex(mod(i,8)+1) cellshape(mod(i,10)+1)];
%             cellind = CellGroup{i};            
%             for j = 1:length(cellind)
%                 t1 = CellTrace_backup(:,cellind(j));
%                 t1 = t1 - mean(t1); 
%                 t = t + t1; %CellTrace_backup(:,cellind(j));
%                 plot(x,t1,ColorIndex(mod(i,8)+1));
%             end
%             plot(x,t/length(cellind),'k')
%             title(num2str(i));
%             % title(q1)
%         end
%          saveas(h,[Ca2Dfile(1:l-4) '_' num2str(le) 'group_Trace.fig'],'fig');
%         
%         
%     end
%     
% end

% 
% qqq = [1 3 4 5 7 9 12];
% figure;
% imagesc(Ca2d); hold on;
% for i = 1:7
%     % subplot(4,5,i)
%     
%     cellind = CellGroup{qqq(i)};
%     q1 = [ColorIndex(mod(i,8)+1) cellshape(mod(i,10)+1)];
%     plot(CellPosition(1,cellind),CellPosition(2,cellind),q1); % ,'linewidth',2);
%     % hold off
%    % title(num2str(i));
% end
%         h = figure;
%         x = 1:Frames;
%         % x = x*0.256;
%         for i = 1:7
%             t = zeros(300,1);            
%             subplot(4,2,i); hold on;
%             q1 = [ColorIndex(mod(i,8)+1) cellshape(mod(i,10)+1)];
%             cellind = CellGroup{qqq(i)};            
%             for j = 1:length(cellind)
%                 t1 = CellTrace_backup(:,cellind(j));
%                 t1 = t1 - mean(t1); 
%                 t = t + t1; %CellTrace_backup(:,cellind(j));
%                 plot(x,t1,ColorIndex(mod(i,8)+1));
%             end
%             plot(x,t/length(cellind),'k')
%             title(num2str(i));
%             % title(q1)
%         end
% 
% CellID = zeros(1,518);
% CellID1 = 1:518;
% st = 1;
% for i = 1: length(CellGroup);
% cellind = CellGroup{i};
% l = length(cellind);
% ed = st+l-1;
% CellID(st:ed) = cellind;
% st = st+l;
% end
% for i = 1:ed
% CellID1(CellID(i)) = 0;
% end
% [e h] =sort(CellID1);
% CellID(ed+1:518) = e(ed+1:518);
% for i = 1:518
% CellTrace(:,i) = CellTrace_backup(:,CellID(i));
% end
% [R P] = corrcoef(CellTrace);
% for i = 1:518
%     R(i,i) = 0;
% end
% figure
% imagesc(R)
% colormap(m)

        
        
        
        
        
        
        
        
        
        
        
% CellID = zeros(1,518);
% CellID1 = 1:518;
% st = 1;
% for i = 1: length(CellGroup);
%     cellind = CellGroup{i};  
%     l = length(cellind);
%     ed = st+l-1;
%     CellID(st:ed) = cellind;
%     st = st+l;
% end
% 
% for i = 1:ed
%     CellID1(CellID(i)) = 0;
% end
%     
% [e h] =sort(CellID1);
% CellID(ed+1:518) = e(ed+1:518);
% 
% for i = 1:518
%     CellTrace(:,i) = CellTrace_backup(:,CellID(i));
% end


%     q1 = [ColorIndex(mod(i,8)+1) cellshape(mod(i,13)+1)];
%     plot(CellPosition(1,cellind),CellPosition(2,cellind),q1,'linewidth',3);
%



% for i = 1:length(CellGroup)
%     q=CellGroup{i};
%     w1 = find(q == 100);
%     w2 = find(q == 198);
%     w3 = find(q == 247);
%     if ~isempty([w1 w2 w3])
%         disp(num2str(i))
%     end
% end
%
%
%
% a1 = [11 41 65 198 200];
% a2 = [66 77 247 100 248];
% figure; imagesc(Ca2d); hold on;
% plot(CellPosition(1,a1),CellPosition(2,a1),'ko','linewidth',3);
% q = CellGroup{20};
% plot(CellPosition(1,q),CellPosition(2,q),'c^','linewidth',2);
%
%
%
%
% CellData.CellPosition = CellPosition;
% CellData.CellArea = CellArea;
% CellData.CellBoundingBox
% % Position of the line scan cells
% LineCellPosition = LineScanData.CellPosition;
% COIPositionX = xst + xLS/2;
%
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