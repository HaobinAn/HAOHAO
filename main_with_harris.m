file_path = 'D:\coner_detection\datasets2\';%change to ur own path
img_path_list=dir(strcat(file_path,'*.JPG'));
img_num=length(img_path_list);
len = size(img_path_list);
for j=1:42;
    
    image_name=img_path_list(j).name;
    image=imread(strcat(file_path,image_name));

    [step,masked]=createMaskfortable(image);%extract color for table

     Nimage = imdilate(step,strel('diamond',8));
     Nimage=segimg(masked,step);%change to binary pic
    figure;
    
    L = bwlabel(Nimage);%mark related area
    stats = regionprops(L);%extract pciture inform
    Ar = cat(1, stats.Area);
    ind = find(Ar ==max(Ar));%find the biggest related area
    Nimage(find(L~=ind))=0;%delete other nosie
%     edgeimg=bwmorph(Nimage,'remove');
    
    Nimage = imresize(Nimage,0.3);%scale it
    
    edgeimg = edge(Nimage,'canny',[0.01,0.2]);

    points=detectHarrisFeatures(edgeimg,'MinQuality',0.3,'FilterSize', 37,'ROI', [618,52,550,500]);
    strongest = selectStrongest(points,2);

    [n,m]=size(edgeimg);
    [H,T,R]=hough(edgeimg,'RhoResolution',0.85,'Theta',-90:1.3:89);
    P  = houghpeaks(H,4);
    lines= houghlines(edgeimg,T,R,P,'FillGap',20,'MinLength',50);%we will have max 4 lines.
    image=imresize(image,0.3);%we scale the orginal pic as the binary pic.
%     imshow(edgeimg);
%     hold on;
    imshow(image);
    hold on;
%     plot(C(:,1),C(:,2),'r-o');
    plot(strongest);
    hold on;
 %%
% rk=zeros(1,4);
% rb=zeros(1,4);
% for a = 1:length(lines)  %plot the lins along the table edges.
%     x1=[lines(a).point1(:,1)];
%     x2=[lines(a).point2(:,1)];
%     y1=[lines(a).point1(:,2)];
%     y2=[lines(a).point2(:,2)];
%     syms k b ;
%     [k,b]=solve(x1*k+b-y1,x2*k+b-y2,k,b);
%     rk(a)=double(k);
%     rb(a)=double(b);
%     x=1:2800;%set the line length  
%     line=x*rk(a)+rb(a);
%     line=int16(line);
%     line1{a}=line;  
% %     plot(line);
% %     hold on;
% end
% %%
% [NL,NW]=size(Nimage);
% for i = 1:length(lines);%solve all the line fuctions
%     k=0;
%    for j=i+1:length(lines);
%          syms xn yn
%          [xn,yn]=solve(xn*rk(i)+rb(i)-yn,xn*rk(j)+rb(j)-yn);%,xn*rk(3)+rb(3)-yn,xn*rk(4)+rb(4)-yn,yn,xn);
%          xnn=double(xn);
%          ynn=double(yn);
%         px(i)=xnn;
%         py(i)=ynn;
%            if 0<xnn && xnn<NW && 0<ynn && ynn<NL;
%                if abs(tan((rk(j)-rk(i))/(1+rk(i)*rk(j))))>tan(10);%夹角太小的点去掉，说明两线重合。 
%                   plot(xnn,ynn,'r-o');
%                   hold on;    
%                end
%            end
%    end
% end
% 
% 
% % pervious=k;
% ppx=px;
% ppy=py;
saveas(gcf,image_name ,'jpg');
close;
end



