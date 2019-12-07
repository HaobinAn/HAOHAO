file_path = 'D:\coner_detection\datasets2\';%change to ur own path
img_path_list=dir(strcat(file_path,'*.JPG'));
img_num=length(img_path_list);
len = size(img_path_list);
for j=21:22;
    image_name=img_path_list(j).name;
    image=imread(strcat(file_path,image_name));
    
    [step,masked]=createMaskfortable(image);%extract color for table
%     [step2,masked2]=createMaskforbook(image);%extract color for table
%     step=step+step2;
%     masked=masked+masked2;
    Nimage=segimg(masked,step);%change to binary pic
    L = bwlabel(Nimage);%mark related area
    stats = regionprops(L);%extract pciture inform
    Ar = cat(1, stats.Area);
    ind = find(Ar ==max(Ar));%find the biggest related area
    Nimage(find(L~=ind))=0;%delete other nosie
%     edgeimg=bwmorph(Nimage,'remove');
     Nimage = imresize(Nimage,0.3);%scale it 
    edgeimg = edge(Nimage,'canny',[0.01,0.2]);
    [n,m]=size(edgeimg);
    [H,T,R]=hough(edgeimg,'RhoResolution',0.85,'Theta',-90:1:89);
    P  = houghpeaks(H,3);%,'threshold',ceil(0.2*max(H(:))));
    lines= houghlines(edgeimg,T,R,P,'FillGap',20,'MinLength',50);%we will have max 4 lines.
    image=imresize(image,0.3);%we scale the orginal pic as the binary pic.
imshow(image);
hold on;
rk=zeros(1,4);
rb=zeros(1,4);
for a = 1:length(lines)  %plot the lins along the table edges.
    x1=[lines(a).point1(:,1)];
    x2=[lines(a).point2(:,1)];
    y1=[lines(a).point1(:,2)];
    y2=[lines(a).point2(:,2)];
    syms k b ;
    [k,b]=solve(x1*k+b-y1,x2*k+b-y2,k,b);
    rk(a)=double(k);
    rb(a)=double(b);
    x=1:2800;%set the line length  
    line=x*rk(a)+rb(a);
    line=int16(line);
    line1{a}=line;  
    plot(line);
    hold on;
end

for i = 1:length(lines);%solve all the line fuctions
   for j=i+1:length(lines);
         syms xn yn
         [xn,yn]=solve(xn*rk(i)+rb(i)-yn,xn*rk(j)+rb(j)-yn);%,xn*rk(3)+rb(3)-yn,xn*rk(4)+rb(4)-yn,yn,xn);
         xnn=double(xn);
         ynn=double(yn);
%          r = 5;
% %         [xx,yy] =rectangle('Position',[xnn-r,ynn-r,2*r,2*r],'Curvature', [1 1]);
%          n=25;
%          xy=zeros(n,2);
%          h=2*pi/(n-1);
%          for q=2:n
%              xy(q,1)=xnn+r*cos((q-1)*h);
%              xy(q,2)=ynn+r*sin((q-1)*h);
% %              sum=L(xy(q-1,1),xy(q-1,2))+L(xy(q,1),xy(q,2));
%                plot(xy(q,1),xy(q,2));
%                hold on;
%          end
         
%          if sum>0;
            plot(xnn,ynn,'r-o');
            hold on;
%          end
   end
end

end



