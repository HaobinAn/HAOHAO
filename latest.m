%change those
original_video = "0.avi";
output_folder = "frames";

video_explode(original_video, output_folder);
a = dir([output_folder '/*.jpg']);
len = size(a,1);

len = len-mod(len,10); %?
img_path_list = dir(strcat((output_folder, ' *.jpg'));

% obj = VideoReader(file_path);
% len = obj.NumberOfFrames;
% 
% len = len-mod(len,10);

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
    image_name=img_path_list(j).name;
    image=imread(strcat(output_folder,image_name));
    %image = read(obj,j);
    
    [step,masked]=createMaskfortable3(image);%extract color for table
    Nimage=segimg(masked,step);%change to binary pic
    
    L = bwlabel(Nimage); %mark related area
    stats = regionprops(L); %extract picture inform
    Ar = cat(1, stats.Area);
    ind = find(Ar ==max(Ar)); %find the biggest related area
    Nimage(find(L~=ind))=0; %delete other noise

    %edges
    edgeimg = edge(Nimage,'approxcanny');
    edgeimglist{picnum}=edgeimg; %for next loop.
    picnum=picnum + 1;

    %find (5) strongest hough lines
    [n,m]=size(edgeimg);
    [H,T,R]=hough(edgeimg,'RhoResolution',1.5,'Theta',-90:1.3:89);
    P  = houghpeaks(H,5);
    lines = houghlines(edgeimg,T,R,P,'FillGap',60,'MinLength',50);

%   figure;
%   imshow(edgeimg);
%   hold on;
%%
    rk = zeros(1,5);
    rb = zeros(1,5);
         
    %find (k, b) for all (max. (5)) lines:         [rk]; [rb]
    for a = 1: length(lines)
            %get the line points
            x1=[lines(a).point1(:,1)];
            x2=[lines(a).point2(:,1)];
            y1=[lines(a).point1(:,2)];
            y2=[lines(a).point2(:,2)];
            
            %calculate (k, b) for every line
            rk(a) = (y2-y1)/(x2-x1);
            syms b;
            if rk(a)==Inf
                rk(a)= 9999999;
            end
                rb(a) = solve(y1 == rk(a)*x1 + b, b);
                
%             %plot the lines (optionally)  
%             plot([x1,x2],[y1,y2], "LineWidth", 2);
%             hold on;
         end
%%
    [NL,NW]=size(Nimage); %len, width of a picture
%%
    %solve all the line fuctions; save intersections:       [cornerarea]
    for ifnf = 1:length(lines)-1
       for jfnf=ifnf+1:length(lines)
             syms xn yn
             [xn,yn]=solve(xn*rk(ifnf)+rb(ifnf)-yn, xn*rk(jfnf)+rb(jfnf)-yn);
             xnn=double(xn);
             ynn=double(yn); 
             if 0<xnn && xnn<NW && 0<ynn && ynn<NL %changed & to && (?) %save intersections that are inside of the picture
                     if abs((rk(jfnf)-rk(ifnf))/(1+rk(ifnf)*rk(jfnf)))>tand(15)  %remove overlapped lines 
                         cornerarea{r}=[xnn,ynn];%save all intersections
                         r=r+1;
                     end
             end
        end
    end
end
if exist('cornerarea','var')
kx=[];
ky=[];
cornersize=size(cornerarea);
for ifc = 1:cornersize(1,2)
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
for ifa=1:length(kx)-1
    for ifb=ifa+1:length(kx)
        if abs(kx(ifa)-kx(ifb))<100 && abs(ky(ifa)-ky(ifb))<100 %remove the points which are too close
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

% % %?????????????ROI
perviousxy{ifxy}=[kkx,kky];
for afr=1:length(kkx);%??ROI?????
    sizex=200;
    sizey=200;
    xback=100;
    yback=100;
    if kkx(afr)-xback+sizex>NW;  
        ROI{afr}=[kkx(afr)-xback,kky(afr)-yback,NW-kkx(afr)+xback-3,sizey];
    elseif kky(afr)-yback+sizey>NL;
        ROI{afr}=[kkx(afr)-xback,kky(afr)-yback,sizex,NL-kky(afr)+yback-3];
    elseif kkx(afr)-xback<0
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
        for b=1:sizeROI(1,2)
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
        saveas(gcf,[image_name,'.jpg']);

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

