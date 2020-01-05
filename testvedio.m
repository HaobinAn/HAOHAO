file_path = 'D:\matlab_workspace\vedio_test\IMG_1753[1].avi';%change to ur own path
obj = VideoReader(file_path);
len = obj.NumberOfFrames;% 读取视频的帧数CurrentTime

len = len-mod(len,10);

r=1;
d=1;
picnum=1;
ifxy=1;
for i=1:10:len
clearvars -except i perviousxy ifxy obj test len;
r=1;
d=1;
picnum=1;
for j=i:i+10
% for j=16:20;
    image = read(obj,j);
    
    [step,masked]=createMaskfortable3(image);%extract color for table
    Nimage=segimg(masked,step);%change to binary pic
    
    L = bwlabel(Nimage);%mark related area
    stats = regionprops(L);%extract pciture inform
    Ar = cat(1, stats.Area);
    ind = find(Ar ==max(Ar));%find the biggest related area
    Nimage(find(L~=ind))=0;%delete other nosie
%     edgeimg=bwmorph(Nimage,'remove');
    
%     Nimage = imresize(Nimage,0.3);%scale it
    
    edgeimg = edge(Nimage,'approxcanny');%,[0.01,0.2]);
    edgeimglist{picnum}=edgeimg;%for next loop.
    picnum=picnum+1;

    [n,m]=size(edgeimg);
    [H,T,R]=hough(edgeimg,'RhoResolution',1.5,'Theta',-90:1.3:89);
    P  = houghpeaks(H,5);
    lines1= houghlines(edgeimg,T,R,P,'FillGap',60,'MinLength',50);%we will have max 4 lines.
%     image=imresize(image,0.3);%we scale the orginal pic as the binary pic.
%     figure;
%     imshow(edgeimg);
%     hold on;
 %%
rk=zeros(1,4);
rb=zeros(1,4);
        for a = 1:length(lines1)  %plot the lins along the table edges.
            x1=[lines1(a).point1(:,1)];
            x2=[lines1(a).point2(:,1)];
            y1=[lines1(a).point1(:,2)];
            y2=[lines1(a).point2(:,2)];
             it(a)=lines1(a).theta;
        end

         lines=lines1;
         for a=1:length(lines);
            x1=[lines(a).point1(:,1)];
            x2=[lines(a).point2(:,1)];
            y1=[lines(a).point1(:,2)];
            y2=[lines(a).point2(:,2)];
            syms k b ;
            [k,b]=solve(x1*k+b-y1,x2*k+b-y2,k,b);
            rk(a)=double(k);
            rb(a)=double(b);
            x=1:2800;
            line=x*rk(a)+rb(a);
            line=int16(line);
%             plot(line);
%             hold on;
         end
%%
[NL,NW]=size(Nimage);
if length(lines)>1;
for ifnf = 1:length(lines)-1;%solve all the line fuctions
   for jfnf=ifnf+1:length(lines);
         syms xn yn
         [xn,yn]=solve(xn*rk(ifnf)+rb(ifnf)-yn,xn*rk(jfnf)+rb(jfnf)-yn);%,xn*rk(3)+rb(3)-yn,xn*rk(4)+rb(4)-yn,yn,xn);
         xnn=double(xn);
         ynn=double(yn); 
         if 0<xnn & xnn<NW & 0<ynn & ynn<NL;
                 if abs((rk(jfnf)-rk(ifnf))/(1+rk(ifnf)*rk(jfnf)))>tand(15);%remove overlapped lines 。 
                  cornerarea{r}=[xnn,ynn];%save all intersections
                  r=r+1;
                 end
         end
    end
end
end
end
if exist('cornerarea','var');
kx=[];
ky=[];
cornersize=size(cornerarea);
for ifc = 1:cornersize(1,2);
%      for j=i+1:cornersize(1,2);
       kx(ifc)=cornerarea{ifc}(1,1);
       ky(ifc)=cornerarea{ifc}(1,2);
%        kx(ifc)=x11;
%        ky(ifc)=y11;
end

%%

% kx=sort(kx);
% ky=sort(ky);
for ifkxa=1:length(kx) 
    kx(kx==0)=[];
end
for ifkya=1:length(kx) 
    ky(ky==0)=[];
end

%%
%使用哈里斯探测指定区域角点
for ifa=1:length(kx)-1
    for ifb=ifa+1:length(kx);
        if abs(kx(ifa)-kx(ifb))<100 && abs(ky(ifa)-ky(ifb))<100;%remove the points which are too closed
            kx(ifa)=0;
            ky(ifa)=0;
        end
    end
end
kkx=nonzeros(kx);
kky=nonzeros(ky);
% 

%     figure;
%     imshow(edgeimg);
%     hold on;
% for ifp=1:length(kkx)
%     plot(kkx(ifp),kky(ifp),'r-o');
%     hold on;
% end

% % %去除相似交点，基于交点设定ROI
perviousxy{ifxy}=[kkx,kky];
for afr=1:length(kkx);%限定ROI不超出图像
    sizex=200;
    sizey=200;
    xback=100;
    yback=100;
    if kkx(afr)-xback+sizex>NW;  
        ROI{afr}=[kkx(afr)-xback,kky(afr)-yback,NW-kkx(afr)+xback-3,sizey];
    elseif kky(afr)-yback+sizey>NL;
        ROI{afr}=[kkx(afr)-xback,kky(afr)-yback,sizex,NL-kky(afr)+yback-3];
    elseif kkx(afr)-xback<0;
        ROI{afr}=[1,kky(afr)-yback,sizex,sizey];
    elseif kky(afr)-yback<0;
        ROI{afr}=[kkx(afr)-xback,1,sizex,sizey];
    elseif kky(afr)-yback<0 && kkx(afr)-xback<0;
        ROI{afr}=[1,1,sizex,sizey];
    elseif kky(afr)-yback+sizey>NL && kkx(afr)-xback+sizex>NW;
        ROI{afr}=[kkx(afr)-xback,kky(afr)-yback,NW-kkx(afr)+xback-3,NL-kky(afr)+yback-3];
    else
        ROI{afr}=[kkx(afr)-xback,kky(afr)-yback,sizex,sizey];
    end
end
picnum=1;   
sizeROI=size(ROI);
for nj=i:i+9
     figure;
        image = read(obj,nj);
%         image=imresize(image,0.3);
%         figure;
%         imshow(edgeimglist{picnum});
        imshow(image);
        hold on;
        for b=1:sizeROI(1,2);
        points=detectHarrisFeatures(edgeimglist{picnum},'MinQuality',0.15,'FilterSize', 47,'ROI',ROI{b});%[1093.18758502544,328.925832619539,70,70]);% ROI{b});
        
        strongest= selectStrongest(points,1);    
        pointnew{d,b}=strongest.Location;
        d=d+1;

        plot(strongest);
        hold on;
        % hold on;
%         clearvars strongest;
        end
        picnum=picnum+1;
        image_name=num2str(nj,'%03d');
        saveas(gcf,['D:\matlab_workspace\vedio_test\image\',image_name,'.jpg']);

        close;
end
else
     for nj=i:i+9
           figure;
           imshow(image);
           hold on;
           image_name=num2str(nj,'%03d');
           saveas(gcf,['D:\matlab_workspace\vedio_test\image\',image_name,'.jpg']);
           close;
     end
end
end




%%
% save as a video
% myObj = VideoWriter('newfile.avi');%初始化一个avi文件
% myObj.FrameRate = 30;
% myObj.Quality=100;
% open(myObj);
% output_path='D:\matlab_workspace\vedio_test\';
% out_path_list=dir(strcat(output_path,'*.JPG'));
% leno=length(out_path_list);
% for frami=1:leno%图像序列个数
%     
%      image_name=out_path_list(frami).name;
%      image=imread(strcat(output_path,image_name));
% %     image=imresize(image,3/10);
% %    image_name=img_path_list(frami).name;
% %    fname=dir(strcat(file_path,'*.JPG'));
% %    frame = imread(fname);
%    writeVideo(myObj,image);
% end
% close(myObj);


