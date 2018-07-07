% setting colormap
clrs=colorcube(12)
clrs(12,:)= [0.9100 0.4100 0.1700] % changing white to orange

%% image sections figure
h1=figure;set(h1,'position',[638 478 1182 172]);
imagesc(mean(spike,3));
colormap jet; minC=0; maxC=1200;  
set(gca,'clim',[minC,maxC]);
for j=1:20
c=xx.(sprintf('s_%d', j));
r=yy.(sprintf('s_%d', j));
hold on;plot(c,r,'w','linewidth',2);
end
axis off
tightfig
%%
h2=figure;set(h2,'position',[638 478 1182 172]);
i=1
for j=1:4   
subplot(1,4,j)
imagesc(mean(subR.(sprintf('s%d', j)),3))
axis off; colormap jet; 
minC=0;maxC=1200;set(gca,'clim',[minC,maxC]);
R=subR.(sprintf('s%d', j));
SectionSize=size(R,2)/3;
mEnd=size(R,2)- SectionSize;
xend=size(R,1);
     for m=0:SectionSize:mEnd  
      rectangle('Position',[m 0.1 SectionSize-2 xend],'EdgeColor',clrs(i,:),'linewidth',2)
    % rectangle('Position',[m 1 SectionSize 31],'EdgeColor','w')
      i=i+1;
     end
end
tightfig
%%
h=figure;% set(h,'position',[810 50 740 736]);
for i=1:size(eo,2)
    subplot(size(eo,2),1,i); plot(bk,eo(:,i),'color',clrs(i,:));axis([-inf inf -inf inf])
    axis([-inf inf -inf inf])
    set(gca,'Visible','off');
end
    ax1=gca;ax1.XAxis.Visible = 'on';



%% 1. FILTER 
% % testing
% 1/(0.002*32)~ 16 frames/sec
%1/(.002*128)=3.90625 ~ 4Hz
% fftshow(eo(:,2),Fs)
%  suptitle( '.0624 HZ')
bk=(1:length(eo))*0.002*128; %or 32 for 16hz and 128 for 4Hz

clear lfp lfp0 lfp1 lfp16 lfp32
%filter at 1/15.625 = 0.0624 HZ cheby1 LPF %FOR 16HZ DATA
%filter at 1/3.90625 = 0.256 HZ cheby1 LPF % FOR4 HZ DATA

for i=1:size(eo,2)
ff=eo(:,i);
[cols k] = cheby1(2,0.256,.15/1,'low');%CHANGE FILTER FOR 4HZ & 16HZ
q = filtfilt(cols,k,ff);
[f1 fftdata1]=fftshow(q,Fs);close
lfp0(:,i)=q;
end
size(lfp0)
%
h=figure;% set(h,'position',[810 50 740 736]);
for i=1:size(eo,2)
    subplot(size(eo,2),1,i); plot(bk,eo(:,i));axis([-inf inf -inf inf])
    title(i,'FontSize',9);
    hold on;p=plot(bk,lfp0(:,i));
    axis([-inf inf -inf inf])
    set(gca,'Visible','off');
end
    ax1=gca;ax1.XAxis.Visible = 'on';
suptitle('lfp filtered 0.0624hZ')
size(lfp0)

%% %% SMOOTH
lfp1=lfp0;

for i=1:size(lfp1,2)
lfp16(:,i)=smooth(lfp1(:,i),16); % 16 datapts=1sec window 
lfp32(:,i)=smooth(lfp1(:,i),32); % 32 datapts=2sec window 
end
% figure;
% bk=(1:length(lfp1(st(j):ed(j),:)))*0.002*32
% subplot(311);plot(bk,lfp1(st(j):ed(j),:)); title('post filtering');axis([-inf inf -inf inf]);
% subplot(312);plot(bk,lfp16(st(j):ed(j),:)); title('soomth 16 datapts= 1 sec'); axis([-inf inf -inf inf]);
% subplot(313);plot(bk,lfp32(st(j):ed(j),:)); title('smooth 32 datapts= 2 sec'); axis([-inf inf -inf inf]);
% % % % % 
% % % % % figure;
% % % % % bk=(1:length(lfp1))*0.002*32
% % % % % subplot(311);plot(bk,lfp1);axis([-inf inf -inf inf]);
% % % % % subplot(312);plot(bk,lfp16); title('soomth 16 datapts= 1 sec'); axis([-inf inf -inf inf]);
% % % % % subplot(313);plot(bk,lfp32); title('smooth 32 datapts= 2 sec'); axis([-inf inf -inf inf]);


% set(h, {'color'}, num2cell(clrs1, 2));
% title('post smoothing ( 1 sec window) ')

lfp=lfp16; size(lfp)
% % % % % % %% plot(lfp) one below the other and together
% % % % % % figure;
% % % % % % lfpFlip=fliplr(lfp);%%NOTE THIS
% % % % % % 
% % % % % % bk=(1:length(lfpFlip))*0.002*32;
% % % % % % p=plot(bk,lfpFlip + 2*repmat(0:(size(lfpFlip,2)-1),size(lfpFlip,1),1),'LineWidth',2);
% % % % % % xlabel('Time (s)');axis([-inf inf -inf inf])
% ax1=gca; set(gca,'Visible','off'); ax1.XAxis.Visible = 'on';    
%%
figure;
for i=1:size(lfp,2)
subplot(12,1,i)    
hold on;plot(bk,lfp(:,i),'LineWidth',1,'color',clrs(i,:));
axis([ -inf inf -inf inf]);%xlabel('sec')
ax1=gca; set(gca,'Visible','off');
% axis off
end
% ax1=gca; set(gca,'Visible','off'); ax1.XAxis.Visible = 'on';    
subplot(12,1,i)    
hold on;plot(bk,lfp(:,i),'LineWidth',2,'color',clrs(i,:));
axis([ -inf inf -inf inf]);%xlabel('sec')
ax1=gca; set(gca,'Visible','on');
set(gca,'Color','none','ytick',[])
% % % % % % 
% % % % % % %% select spikes size
% % % % % % fps=16;
% % % % % % [pk,l,w,p]=findpeaks(b,fps,'MinPeakHeight',40,'MinPeakDistance',20)
% % % % % % figure;findpeaks(b,fps,'MinPeakHeight',40,'MinPeakDistance',20,'Annotate','extents')
% % % % % % 
% % % % % % 
% % % % % % size(pk)
% % % % % % 
% % % % % % %% spike st and ed  
% % % % % % cr=.002*32;% for 16fs: 32 pixels/line and 0.02  ?
% % % % % % st=(l-15)/cr 
% % % % % % ed=(l+25)/cr
% % % % % % 
% % % % % % %% sanity check st and ed
% % % % % % figure;plot(lfp(:,1:3));
% % % % % % hold on; plot(lfp(:,6))
% % % % % % for i=1:6
% % % % % % h=vline(st(i),'m')
% % % % % % h=vline(ed(i),'k')
% % % % % % end
% % % % % % 
% % % % % % %% fix if needed
% % % % % %  [x y]=ginput(4)
% % % % % % st(2)=x(1)
% % % % % % ed(2)=x(2)
% % % % % % st(3)=x(3)
% % % % % % ed(3)=x(4)
% % % % % % 
% % % % % % %% FINDpeaks and FWHM as rise time
% % % % % % % clear rtL rtV ptL ptV r2L p2L r2V p2V  %ephysST
% % % % % % 
% % % % % % StPk=2;%SELECT starting pk number-> SPIKE NUMBERS TO PROCESS!!
% % % % % % 
% % % % % % for i=2:11 % REGION#
% % % % % %     figure; subplot(4,4,[1,2]);plot(lfp(:,i));
% % % % % %     axis([-inf inf -inf inf]);title(['Region ' num2str(i)])       
% % % % % %     
% % % % % %     for j=StPk:size(pk,1) %SPIKE#
% % % % % %         % SET  region 1 nad 12 values to zero!!!
% % % % % %        sp=lfp(st(j):ed(j),i);
% % % % % %        [spk,sl,sw,spr] =findpeaks(sp,'MinPeakHeight',.8,'MinPeakWidth',5,'Annotate','extents');
% % % % % %        [a1 a2]=max(spk);
% % % % % %         if  ~isempty(a2)
% % % % % %             ptV(:,i,j)=a1;
% % % % % %             ptL(:,i,j)=sl(a2);
% % % % % %             rtV(:,i,j)=max(spk)-spr(a2)/2;
% % % % % %             [minDiffV, indMin] = min(abs(sp((ptL(:,i,j)-(sw(a2)):ptL(:,i,j)))- rtV(:,i,j)));
% % % % % %             rtL(:,i,j)=indMin+(ptL(:,i,j)-(sw(a2)))-1 ; 
% % % % % %             subplot(4,4,j+2)
% % % % % %             bk=(1:length(sp))*0.002*32;
% % % % % %             plot(bk,sp);title(['spike ' num2str(j)]);axis([-inf inf -inf inf])       
% % % % % %                 hold on
% % % % % %             line([rtL(:,i,j)*cr rtL(:,i,j)*cr],[min(sp)  max(sp)],'color','k'); %h=vline(rtL(:,i,j)*cr,'m');
% % % % % %             line([ptL(:,i,j)*cr ptL(:,i,j)*cr],[min(sp) max(sp)],'color','g'); % h=vline(ptL(:,i,j)*cr,'k');
% % % % % %        else
% % % % % %           rtL(:,i,j)=0;
% % % % % %           ptL(:,i,j)=0;
% % % % % %           ptV(:,i,j)=0;
% % % % % %           rtV(:,i,j)=0;
% % % % % %           subplot(4,4,j+2);
% % % % % %           bk=(1:length(sp))*0.002*32;
% % % % % %           plot(bk,sp);title(['spike ' num2str(j)])       
% % % % % %           axis([-inf inf -inf inf])       
% % % % % %         end
% % % % % %     end
% % % % % % end
% % % % % % 
% % % % % % %% % SET  region 1 nad 12 values to zero!!!
% % % % % % rtL(:,12,j)=0;
% % % % % % ptL(:,12,j)=0;
% % % % % % ptV(:,12,j)=0;
% % % % % % rtV(:,12,j)=0;
% % % % % % 
% % % % % % 
% % % % % % % rise time nad peak time in seconds
% % % % % % rtL1=round(squeeze(rtL)*cr,2)
% % % % % % ptL1=round(squeeze(ptL)*cr,2)
% % % % % % 
% % % % % % %% That figure in paper?
% % % % % % c=clrs;
% % % % % % % figure
% % % % % % h=figure(42)
% % % % % % subplot(311)
% % % % % % clear rtALL ptALL
% % % % % % rtALL=rtL1;
% % % % % % rtALL(rtALL==0)= nan;
% % % % % % 
% % % % % % for i=1:size(rtALL,2)
% % % % % % x=rtALL(:,i)'
% % % % % % y=1:12
% % % % % % I = ~isnan(x) & ~isnan(y);hold on
% % % % % % plot(x(I),y(I),'k')
% % % % % % end
% % % % % % 
% % % % % % for i=1:size(rtALL,2)
% % % % % % hold on; 
% % % % % % scatter(rtALL(:,i),y,35,c, 'filled')
% % % % % %    
% % % % % % end
% % % % % % 
% % % % % % locx=7;%% SELECT
% % % % % % for i=[2 6]%1:size(rtALL,2)
% % % % % % x1 =(rtALL(locx,i)+.1)
% % % % % % y1 =locx;
% % % % % % txt = [num2str(i)]
% % % % % % text(x1,y1,txt,'Color','k','FontSize',10)
% % % % % % end
% % % % % % % view([90 -90]) %// instead of normal view, which is view([0 90])
% % % % % % xlabel('FWHM: rise time (s)');ylabel('Region ')
% % % % % % % h=hline(9,'k');h=hline(6,'k');h=hline(3,'k');
% % % % % % grid on; ax=gca; ax.YTick=([1:12]);
% % % % % % axis([-inf inf 1 12]);title('spike#')
% % % % % % % grid minor
% % % % % % 
% % % % % % 
% % % % % % %% that figure for peaks
% % % % % % 
% % % % % % h=figure(42)
% % % % % % subplot(313)
% % % % % % ptALL=ptL1;
% % % % % % ptALL(ptALL==0)= nan;
% % % % % % 
% % % % % % for i=1:size(ptALL,2)
% % % % % % x=ptALL(:,i)'
% % % % % % y=1:12
% % % % % % I = ~isnan(x) & ~isnan(y);hold on
% % % % % % plot(x(I),y(I),'k')
% % % % % % end
% % % % % % 
% % % % % % c=clrs;
% % % % % % for i=1:size(ptALL,2)
% % % % % % hold on; 
% % % % % % scatter(ptALL(:,i),y,35,c, 'filled')
% % % % % % end
% % % % % % 
% % % % % % % locx=10;%% SELECT
% % % % % % for i=1:size(ptALL,2)
% % % % % % x1 =(ptALL(locx,i))+.1
% % % % % % y1 =locx;
% % % % % % txt = [num2str(i)]
% % % % % % text(x1,y1,txt,'Color','r','FontSize',10)
% % % % % % end
% % % % % % % view([90 -90]) %// instead of normal view, which is view([0 90])
% % % % % % xlabel('peak time (s)');ylabel('Region ')
% % % % % % % h=hline(9,'k');h=hline(6,'k');h=hline(3,'k');
% % % % % % grid on; ax=gca; ax.YTick=([1:12]);
% % % % % % axis([-inf inf 1 12]);title('spike#')
% % % % % % % grid minor
% % % % % % 
% % % % % %  
% % % % % % 
% % % % % % 
% % % % % % 
% % % % % % %% figure out Ephys psike starting- manualyy
% % % % % % % clear ephysST
% % % % % % Fs=10000; 
% % % % % % j=6
% % % % % %        eR=(ephys1((st(j))*cr*Fs:ed(j)*cr*Fs));
% % % % % %         eF=(q2((st(j))*cr*Fs:ed(j)*cr*Fs));
% % % % % %         fi=(1:length(eR))/(cr*Fs);max(fi)
% % % % % %         figure; plot(fi,eR);hold on;plot(fi,eR);
% % % % % %         plot(fi,eF*2,'color','k');
% % % % % %        sp=lfp((st(j)):ed(j),1:3);
% % % % % %        plot(sp/5,'g')
% % % % % %        sp=lfp((st(j)):ed(j),7);
% % % % % %        plot(sp/5,'m')
% % % % % %        
% % % % % %         [x y]=ginput(1)
% % % % % %         ephysST(:,j)=x;
% % % % % %         h=vline( ephysST(:,j),'m');
% % % % % %         pl=round(mean(ptL(:,:,j)))
% % % % % %         h=vline(ptL(:,i,j),'r');
% % % % % % 
% % % % % % %% 2*RMS ( average 3)
% % % % % % 
% % % % % % for i=1:size(lfp,2) % REGION
% % % % % %             sample=lfp(:,i);
% % % % % %             h=figure;%set(h,'position',[790 70 968 461]);
% % % % % %             subplot(3,4,[1,2,3,4])
% % % % % %             plot(sample); hold on; title([' region=' num2str(i)])
% % % % % %             Z = smooth(sample,55);plot(Z);hold on
% % % % % % %             ZT=detrend(Z);
% % % % % % %             plot(ZT)
% % % % % % %             Z=ZT
% % % % % %             h=hline(min(sample),'k');
% % % % % %     if i<4 | i>10
% % % % % %                 axis([-inf inf -0.05 .3])
% % % % % %     else
% % % % % %     end            
% % % % % % %           MANUAL
% % % % % %             [x y]=ginput(6);cw(:,i)=x;
% % % % % % %           WITH SAVED DATA
% % % % % %             x=cw(:,i);
% % % % % %             a1=mean(Z((x(1):x(2))));a2=mean(Z((x(3):x(4))));
% % % % % %             a3=mean(Z((x(5):x(6))));
% % % % % %             mu=mean([a1 a2 a3]);
% % % % % %             RMS = rms(mu);2*RMS
% % % % % %             for ii=1:6
% % % % % %             h= vline(x(ii),'k');
% % % % % %             end
% % % % % % %             xlim([-inf inf])
% % % % % %            
% % % % % % 
% % % % % %     for j=StPk:size(pk,1)%SELECT pk-> SPIKE NUMBERS TO PROCESS!!
% % % % % % %  j=3;
% % % % % % 
% % % % % %         sp=lfp(st(j):ed(j),i);
% % % % % %         subplot(3,4,j+4)
% % % % % %         plot(sp);hold on;title([' spike=' num2str(j)]);
% % % % % %         Z = smooth(sp,25);plot(Z);
% % % % % %         h=vline(rtL(:,i,j),'k');
% % % % % %         h=vline(ephysST(:,j),'g');
% % % % % %         
% % % % % %         [spk,sl,sw,spr] =findpeaks(sp,'MinPeakHeight',.06,'MinPeakWidth',5,'Annotate','extents');
% % % % % % %         figure;findpeaks(sp,'MinPeakHeight',.06,'MinPeakWidth',5,'Annotate','extents');
% % % % % %         [a1 a2]=max(spk);
% % % % % % 
% % % % % %     if ~isempty(a2)
% % % % % %         hold on; plot(sl(a2),spk(a2),'o');
% % % % % %         p2L(:,i,j)=sl(a2); 
% % % % % %         p2V(:,i,j)=spk(a2);
% % % % % %         r2V(:,i,j)=2*RMS;
% % % % % %         ee=sl(a2); % peak location
% % % % % %             if ephysST(:,j)>(ee-sw(a2)/2+10)
% % % % % %             window=sl(a2)-ephysST(:,j)
% % % % % %             else
% % % % % %             window=sw(a2)/2+10
% % % % % %             end
% % % % % %                 ss=ee-window;
% % % % % %                 ss(ss<0)=1;
% % % % % %                 [md im]=min(abs(round(sp(ss:ee)-(2*RMS),4)));
% % % % % %                 im+ss-1 
% % % % % %                 r2L(:,i,j)=im+ss-1 ; 
% % % % % %                 hold on; h=hline(mu,'k');
% % % % % %                 h=hline(2*RMS,'k');
% % % % % %                 h=vline(r2L(:,i,j),'r'); %%RISE TIME- loc of 2* RMS
% % % % % %                 axis([-inf size(sp,1) min(sp) max(sp)]);
% % % % % %     else
% % % % % %         r2L(:,i,j)=0;
% % % % % %         r2V(:,i,j)=0;
% % % % % %         p2L(:,i,j)=0;
% % % % % %         p2V(:,i,j)=0;
% % % % % %     end
% % % % % %     
% % % % % %     end
% % % % % % 
% % % % % % end
% % % % % % 
% % % % % % r2L1=round(squeeze(r2L)*cr,2)
% % % % % % 
% % % % % % 
% % % % % % %% THAT FIGURE RMS
% % % % % % h=figure(42)
% % % % % % subplot(312)
% % % % % % r2ALL=r2L1;
% % % % % % r2ALL(r2ALL==0)= nan;
% % % % % % 
% % % % % % for i=1:size(r2ALL,2)
% % % % % % x=r2ALL(:,i)'
% % % % % % y=1:12
% % % % % % I = ~isnan(x) & ~isnan(y);hold on
% % % % % % plot(x(I),y(I),'k')
% % % % % % end
% % % % % % 
% % % % % % c=clrs;
% % % % % % for i=1:size(r2ALL,2)
% % % % % % hold on; 
% % % % % % scatter(r2ALL(:,i),y,35,c, 'filled')
% % % % % % end
% % % % % % 
% % % % % % % locx=10;%% SELECT
% % % % % % for i=1:size(r2ALL,2)
% % % % % % x1 =(r2ALL(locx,i)+.1)
% % % % % % y1 =locx;
% % % % % % txt = [num2str(i)]
% % % % % % text(x1,y1,txt,'Color','k','FontSize',10)
% % % % % % end
% % % % % % view([90 -90]) %// instead of normal view, which is view([0 90])
% % % % % % xlabel('2*RMS : rise time (s) ');ylabel('Region ')
% % % % % % % h=hline(9,'k');h=hline(6,'k');h=hline(3,'k');
% % % % % % grid on; ax=gca; ax.YTick=([1:12]);
% % % % % % axis([-inf inf 1 12]);title('spike#')
% % % % % % % grid minor
% % % % % % 
% % % % % % %% PLOT  individual spikes+ raw ephys TOGETHER
% % % % % % figure;
% % % % % % st1=st%+180;
% % % % % % ed1=ed%-100;
% % % % % % for j=StPk:size(pk,1) %SELECT pk-> SPIKE NUMBERS TO PROCESS!!
% % % % % %     subplot(5,1,j-1); hold on
% % % % % %     Fs=10000;
% % % % % %     eR=(ephys1(st1(j)*cr*Fs:ed1(j)*cr*Fs));
% % % % % %     fi=(1:length(eR))/(Fs);max(fi)
% % % % % %     eF=(q2(st1(j)*cr*Fs:ed1(j)*cr*Fs));
% % % % % %     fi=(1:length(eF))/(Fs);max(fi)
% % % % % % 
% % % % % %     plot(fi,eR);hold on; 
% % % % % %     plot(fi,eF,'color',[.6 0.3 .9]);hold on;
% % % % % % %     clrs=colorcube(12);
% % % % % %         for i=1:12
% % % % % %             sp=lfp(st1(j):ed1(j),i);
% % % % % %             bk=(1:length(sp))*0.002*32;;
% % % % % %             plot(bk,sp/3,'color',clrs(i,:));hold on
% % % % % %         end
% % % % % %     axis([bk(1)+10 bk(end)-10 -inf inf]);title(['spike# ' num2str(j)],'FontSize',8)
% % % % % % ax1=gca; set(gca,'Visible','off'); ax1.XAxis.Visible = 'on';    
% % % % % % 
% % % % % % end
% % % % % % xlabel('sec');
% % % % % % legend('raw ephys','filterd ephys')
% % % % % % 
% % % % % % %%
% % % % % % 
% % % % % % h=figure
% % % % % % % set(h,'position',[1 524 1525 274]);
% % % % % % 
% % % % % % 
% % % % % % for j=[2 ]
% % % % % %         hold on
% % % % % %         eR=(ephys1(st(j)*cr*Fs:ed(j)*cr*Fs));
% % % % % %         fi=(1:length(eR))/(Fs);max(fi)
% % % % % %         plot(fi,eR*10);hold on  
% % % % % %         ylim([-inf inf])
% % % % % % 
% % % % % %     for i=1:11
% % % % % %         
% % % % % %     bk=(1:length(lfp(st(j):ed(j),i)))*0.002*32;
% % % % % %      plot(bk,lfp(st(j):ed(j),i),'color',clrs(i,:));hold on
% % % % % %      line([r2L(:,i,j)*cr r2L(:,i,j)*cr],[-.2 1.8],'Color',clrs(i,:));
% % % % % % % xlim([ptL(:,i,j)*cr-100*cr ptL(:,i,j)*cr+2*cr])
% % % % % % % xlim([ptL(:,i,j)*cr-100*cr 30])
% % % % % % 
% % % % % % xlabel('sec')
% % % % % % ax1=gca; set(gca,'Visible','off'); ax1.XAxis.Visible = 'on';    
% % % % % %     end
% % % % % % end
% % % % % % 
% % % % % % suptitle('Rise time: 2*RMS')
% % % % % % %%
% % % % % % 
% % % % % % h=figure
% % % % % % set(h,'position',[1 524 1525 274]);
% % % % % % for j=2:6
% % % % % %         subplot(1,6,j);hold on
% % % % % %     for i=4:11
% % % % % %         
% % % % % %        bk=(1:length(lfp(st(j):ed(j),i)))*0.002*32;
% % % % % %      plot(bk,lfp(st(j):ed(j),i),'color',clrs(i,:));
% % % % % %      line([rtL(:,i,j)*cr rtL(:,i,j)*cr],[-.2 1.8],'Color',clrs(i,:));
% % % % % % % line([rtL(:,i,j) rtL(:,i,j)],[-.2 1.8],'Color','g');
% % % % % % xlim([ptL(:,i,j)*cr-100*cr ptL(:,i,j)*cr+2*cr])
% % % % % % ylim([-inf inf])
% % % % % % xlabel('sec')
% % % % % % ax1=gca; set(gca,'Visible','off'); ax1.XAxis.Visible = 'on';    
% % % % % %     end
% % % % % % end
% % % % % % suptitle('Rise time: FWHM')
% % % % % % 
% % % % % % 
% % % % % % %% plot average of all rise times together
% % % % % % 
% % % % % % figure;plot(nanmean(rtALL'),1:12);hold on; plot(nanmean(r2ALL'),1:12)
% % % % % % axis([-inf inf 1 12])
% % % % % % grid on
% % % % % % grid minor
% % % % % % legend('Rise time by FWHM','Riste time by 2*RMS')
% % % % % % ylabel('region#')
% % % % % % title( 'average rise time over all regions ');xlabel('sec')
% % % % % % 
% % % % % % %%
% % % % % % % %% plot individual spikes with rt pt and r2
% % % % % % % %  figure
% % % % % % % %     m=1;   
% % % % % % % for j=StPk:size(pk,1) %SELECT pk-> SPIKE NUMBERS TO PROCESS!!
% % % % % % % 
% % % % % % %     figure
% % % % % % %     m=1;
% % % % % % % %   subplot(4,4,1);plot(lfp);axis([-inf inf -inf inf])
% % % % % % %     subplot(3,6,[1,2,3,4,5,6]);
% % % % % % %     hold on
% % % % % % %     Fs=10000;
% % % % % % %     eR=(ephys1(st(j)*cr*Fs:ed(j)*cr*Fs));
% % % % % % %     fi=(1:length(eR))/(Fs);max(fi)
% % % % % % % %   eF=(q2(st(j)*cr*Fs:ed(j)*cr*Fs));
% % % % % % %     plot(fi,eR);hold on; 
% % % % % % %     for i=1:12
% % % % % % %             sp=lfp(st(j):ed(j),i);
% % % % % % %             bk=(1:length(sp))*0.002*32;;
% % % % % % %             plot(bk,sp,'color',clrs(i,:));hold on
% % % % % % %     end
% % % % % % %    axis([-inf inf -inf inf]);title(['spike# ' num2str(j)])
% % % % % % %    xlabel('sec')
% % % % % % % 
% % % % % % % %    axis off
% % % % % % %     
% % % % % % %     
% % % % % % %     for i=1:size(lfp,2)
% % % % % % %         sp=lfp(st(j):ed(j),i);
% % % % % % % %         % figure;
% % % % % % % %         subplot(size(pk,1),13,m+1)
% % % % % % %         subplot(3,6,m+6)
% % % % % % % 
% % % % % % %         plot(sp,'color',clrs(i,:));hold on
% % % % % % %         title(['R' num2str(i)],'FontSize',7)
% % % % % % % % %         Z = smooth(sp,25);plot(Z)
% % % % % % %         h=vline(rtL(:,i,j),'k');
% % % % % % %         h=vline(r2L(:,i,j),'m');
% % % % % % %         h=vline(ephysST(:,j),'g');
% % % % % % %         m=m+1
% % % % % % % %         axis([200 ed(j)/10+200 min(sp)-.01 max(sp)]);
% % % % % % %     end
% % % % % % %     m=m+1;
% % % % % % %     
% % % % % % % end
% % % % % % 
% % % % % % %%
% % % % % % % %% save
% % % % % % % save rtL1
% % % % % % % save ptL1
% % % % % % % save l
% % % % % % % save lfp
% % % % % % % save xx
% % % % % % % save yy
% % % % % % % save BW1
% % % % % % % 
% % % % % % % %% %% FIND rise time 2*RMS of callbearation window
% % % % % % % % 
% % % % % % % j=3;i=3;
% % % % % % % %selecting manualy a20 to 30 datapts as CW 
% % % % % % % m=1
% % % % % % % for j=1:size(lfp,2)
% % % % % % %     sp=lfp(st(j):ed(j),:);
% % % % % % %     figure;plot(sp)
% % % % % % %    [x,y] = ginput(2);
% % % % % % %    cwST(m)=x(1)
% % % % % % %    cwED(m)=x(2)
% % % % % % %    m=m+1
% % % % % % % end
% % % % % % % %
% % % % % % % % drawing mean and RMS
% % % % % % % j=5
% % % % % % % for i=1:size(lfp,2)
% % % % % % %     sp=lfp(st(j):ed(j),i);
% % % % % % %     figure;plot(sp);title(i)
% % % % % % % 
% % % % % % % h=vline(round(cwED(j)),'k');
% % % % % % % h=vline(round(cwST(j)),'k');
% % % % % % % cw=sp(round(cwST(j)):round(cwED(j)));
% % % % % % % RMS = sqrt(nanmean(cw.^2));
% % % % % % % r2V(i)=2*RMS;  
% % % % % % % % r3=3*RMS;
% % % % % % % mu=nanmean(cw);
% % % % % % % h=hlin1e(mu,'k');
% % % % % % % h=hline(r2V(i),'k');
% % % % % % % % h=hline(r3,'k');
% % % % % % % m=((sp(cwED(j):ptL(:,i,j))- r2V(i)));
% % % % % % % [minDiffV, indMin]=min(m(m>0));
% % % % % % % r2L(i)=cwED(j)+indMin-1;
% % % % % % % h=vline(r2L(i),'g');
% % % % % % % 
% % % % % % % 
% % % % % % %     
% % % % % % % end
% % % % % % %%
% % % % % % 
% % % % % % 
% % % % % % %% plot sequence with time diff as text in image and plots
% % % % % % %    plot that graph of rise times
% % % % % % % 
% % % % % % % figure;hold on
% % % % % % % for j=1:6
% % % % % % %     plot(rtL(:,j,2),'o')
% % % % % % %     view([90 -90])
% % % % % % % end
% % % % % % %   
% % % % % %  