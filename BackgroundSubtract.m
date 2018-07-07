%% C. PREPPING 

% C1. bACKGROUND SUbTRACT
spike=CaData;
h=figure;imagesc(mean(double(spike),3))
set(h,'position',[719   455   831   238]);
colormap jet
minC=0; %round(min(min(df)))  
maxC=500; %round(max(max(df)))  
set(gca,'clim',[minC,maxC]);

%select background
clear bm bm1 bm2 spike_bs
bm=(mean(double(spike),3));
    for m=1:3 
        rows=getrect
        bm1=bm(rows(2):rows(2)+rows(4),rows(1):rows(1)+rows(3));
        % figure;imagesc(bm1)
        bm2(m)=mean(mean(bm1))
    end
background_mean=mean(bm2);

% background subtract
    for m=1:size(spike,3);
        spike_bs(:,:,m)= (spike(:,:,m))-(background_mean);
    end;
% spike_bs(spike_bs<0) = 0;
% 
% close; h=figure; %set(h,'position',[[1000,1161,1013,176]]);
% subplot(211);imagesc(mean(double(spike),3));title('raw')
% subplot(212);imagesc(mean(double(spike_bs),3)) ;title('background subtracted')

%
h=figure;imagesc(mean(double(spike_bs),3))
set(h,'position',[719   455   831   238]);
colormap jet
minC=0; %round(min(min(df)))  
maxC=500; %round(max(max(df)))  
set(gca,'clim',[minC,maxC]);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
 
% figure;  imagesc(squeeze(sum(CaData,3))); hold on;

CellCount = 0;
[c f]=butter(2,5/500,'high');
m=1;
for i = 1:32
    t = smooth(w(i,:),5);
    t2 = filtfilt(c,f,t);
    [PEK LOC] = findpeaks(t2,'MinPeakHeight',1e5, 'MinPeakDistance',20);
    t3 = -t2;
    [PEK1 LOC1] = findpeaks(t3,'MinPeakWidth',3);
    [PEK2 LOC2] = findpeaks(-t,'MinPeakWidth',5);
    LOCV = sort([LOC1' LOC2']);
    if ~isempty(LOC)
%         figure
%         plot(t)
%         hold on
%         plot(LOC, t(LOC),'r*');
%         plot(LOC1,t(LOC1),'bo');
        for j = 1:length(LOC)
            st = LOCV(find(LOCV<LOC(j),1,'Last'));
            if isempty(st)
                st =1;
            end
            
            ed = LOCV(find(LOCV>LOC(j),1,'first'));
            if isempty(ed)
                ed =512;
            end
                
            CellCount = CellCount + 1;
            nt = squeeze(sum(CaData(i,st:ed,:)));
            CellTrace(CellCount,:) = nt';
            CellLocation(CellCount,1:3) = [i st ed];
            hold on;plot([st ed ed st st],[i-0.5 i-0.5 i+0.5 i+0.5 i-0.5],'color',clrs(m,:),'LineWidth',2);m=m+1;
%             plot([st ed ed st st],[i-0.5 i-0.5 i+0.5 i+0.5 i-0.5],'w','LineWidth',2);m=m+1;
m
        end
    end
end

% 
% for k = 1:CellCount
%     st = CellLocation(k,2);
%     ed = CellLocation(k,3);
%     i = CellLocation(k,1);
%     plot([st ed ed st st],[i-0.5 i-0.5 i+0.5 i+0.5 i-0.5],'w');
% end
