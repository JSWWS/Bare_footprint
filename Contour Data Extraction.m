clear all;
clc;
%% import data
List_all =dir('F:\XXX\*.jpg');
len_all =length(List_all);
for iall=1:len_all
    file_name_all{iall}=List_all(iall).name;
    img1=importdata(file_name_all{iall}); 
    %% power-law transformed
    file_name_all_origin=file_name_all{iall};
    file_name_all_change=file_name_all{iall};
    file_name_all_change(end-3:end)=[];
    filename_zhi=['zhi_',file_name_all_change];
    filename_gen=['gen_',file_name_all_change];
    C = 5;
    X = 20;
    f = myImageMiLv(img1,C,X);
    
    %% Adaptive median filtering
    
    image_gray=rgb2gray(f);
    ff =image_gray;
    alreadyProcessed = false(size(image_gray));
    Smax=7;
    for k = 3:2:Smax
        zmin = ordfilt2(image_gray, 1, ones(k, k), 'symmetric');
        zmax = ordfilt2(image_gray, k * k, ones(k, k), 'symmetric');
        zmed = medfilt2(image_gray, [k k], 'symmetric');
        
        processUsingLevelB = (zmed > zmin) & (zmax > zmed) & ~alreadyProcessed;
        zB = (image_gray > zmin) & (zmax > image_gray);
        outputZxy  = processUsingLevelB & zB;
        outputZmed = processUsingLevelB & ~zB;
        ff(outputZxy) = image_gray(outputZxy);
        ff(outputZmed) = zmed(outputZmed);
        
        alreadyProcessed = alreadyProcessed | processUsingLevelB;
        if all(alreadyProcessed(:))
            break;
        end
    end
    ff(~alreadyProcessed) = zmed(~alreadyProcessed);
    
    %% Filling and corrosion
 
    A3(:,:,1)=ff;
    A3(:,:,2)=ff;
    A3(:,:,3)=ff;
    thresh = graythresh(A3);
    BW= im2bw(A3,thresh);

    IBW = ~BW;
    P = max(cell2mat(struct2cell(regionprops(~BW,'area'))));
    BW2 = bwareaopen(IBW,round(P/3),8);
    F1 = imfill(BW2,'holes');
    SE = strel('disk',3);
    F2 = imdilate(F1,SE,'same');
    BW3 = bwperim(F2);
    BW4 = ~BW3;
    [yRows, xColumns] = find(BW3);
    mediany=find(yRows==round(median(yRows)));
    startdot(1,2)=xColumns(mediany(1,1),1);
    startdot(1,1)=yRows(mediany(1,1),1);
    dotorigin_origin = bwtraceboundary(BW3,startdot,'W',8,Inf,'counterclockwise');
    dotorigin_origin(:,[1 2]) = dotorigin_origin(:,[2 1]);
    windowWidth = 45;
    polynomialOrder = 3;
    dotorigin_origin = sgolayfilt(dotorigin_origin, polynomialOrder, windowWidth);
    
    bottomy=find(yRows==max(yRows));
    bottomx=find(xColumns==xColumns(bottomy(1,1),1));
    startdot(1,2)=xColumns(bottomx(length(bottomx)-1,1),1);
    startdot(1,1)=yRows(bottomx(length(bottomx)-1,1),1);
    dotorigin_origin_bottom = bwtraceboundary(BW3,startdot,'W',8,Inf,'counterclockwise');
    crosspoint=2;
    while length(dotorigin_origin_bottom)<1000
        startdot(1,2)=xColumns(bottomx(length(bottomx)-crosspoint,1),1);
        startdot(1,1)=yRows(bottomx(length(bottomx)-crosspoint,1),1);
        dotorigin_origin_bottom = bwtraceboundary(BW3,startdot,'W',8,Inf,'counterclockwise');
        crosspoint=crosspoint+1;
    end
    dotorigin_origin_bottom(:,[1 2]) = dotorigin_origin_bottom(:,[2 1]);
    dotorigin_origin_bottom = sgolayfilt(dotorigin_origin_bottom, polynomialOrder, windowWidth);
    
    
    %% angle bisector
    Convex_hull %Convex hull algorithm
    k1=(plot2(1,2)-plot1(1,2))/(plot2(1,1)-plot1(1,1));
    b1=(-1)*k1*plot1(1,1)+plot1(1,2);
    k2=(plot4(1,2)-plot3(1,2))/(plot4(1,1)-plot3(1,1));
    b2=(-1)*k2*plot3(1,1)+plot3(1,2);
    [Xcross,Ycross]=linecross(k1,b1,k2,b2);
    close all
    
    k1angle=(atan(k1))*(180/pi);
    if k1angle<0
        k1angle=180+k1angle;
    end
    k2angle=(atan(k2))*(180/pi);
    if k2angle<0
        k2angle=180+k2angle;
    end
    angle_divide=(k1angle+k2angle)/2;
    kdivide=tan((angle_divide/180)*pi);
    b_dline=Ycross-kdivide*Xcross;
    x_dline1=(min(dotorigin_origin(:,2))-b_dline)/kdivide; 
    x_dline2=(max(dotorigin_origin(:,2))-b_dline)/kdivide;

    
    %% angle of rotation
    rotateangle = angle_divide-90;
    xrotatecenter=(x_dline2+x_dline1)/2;
    yrotatecenter=(max(dotorigin_origin(:,2))+min(dotorigin_origin(:,2))+100)/2;
    
    rotateAroundd=rotateAround(BW3,yrotatecenter,xrotatecenter,rotateangle);
    for idot=1:length(dotorigin_origin)
        dotorigin_origin_rotate(idot,1) = (dotorigin_origin(idot,1)-xrotatecenter)*cosd(rotateangle) + (dotorigin_origin(idot,2)-yrotatecenter)*sind(rotateangle) + xrotatecenter;
        dotorigin_origin_rotate(idot,2) = -(dotorigin_origin(idot,1)-xrotatecenter)*sind(rotateangle) + (dotorigin_origin(idot,2)-yrotatecenter)*cosd(rotateangle) + yrotatecenter;
    end
    
    %% The coordinates of the point selection after rotation
    
    Xplot1 = plot1(1,1);
    Yplot1 = plot1(1,2);
    angplot1 = rotateangle;
    Xcplot1 = xrotatecenter;
    Ycplot1 = yrotatecenter;
    Xrotateplot1 =  (Xplot1-Xcplot1)*cosd(angplot1) + (Yplot1-Ycplot1)*sind(angplot1) + Xcplot1;
    Yrotateplot1 = -(Xplot1-Xcplot1)*sind(angplot1) + (Yplot1-Ycplot1)*cosd(angplot1) + Ycplot1;
    plot1rotate = [Xrotateplot1,Yrotateplot1];
    
    Xplot2 = plot2(1,1);
    Yplot2 = plot2(1,2);
    angplot2 = rotateangle; 
    Xcplot2 = xrotatecenter;
    Ycplot2 = yrotatecenter;
    Xrotateplot2 =  (Xplot2-Xcplot2)*cosd(angplot2) + (Yplot2-Ycplot2)*sind(angplot2) + Xcplot2;
    Yrotateplot2 = -(Xplot2-Xcplot2)*sind(angplot2) + (Yplot2-Ycplot2)*cosd(angplot2) + Ycplot2;
    plot2rotate = [Xrotateplot2,Yrotateplot2];
    
    Xplot3 = plot3(1,1);
    Yplot3 = plot3(1,2);
    angplot3 = rotateangle;
    Xcplot3 = xrotatecenter;
    Ycplot3 = yrotatecenter;
    Xrotateplot3 =  (Xplot3-Xcplot3)*cosd(angplot3) + (Yplot3-Ycplot3)*sind(angplot3) + Xcplot3;
    Yrotateplot3 = -(Xplot3-Xcplot3)*sind(angplot3) + (Yplot3-Ycplot3)*cosd(angplot3) + Ycplot3;
    plot3rotate = [Xrotateplot3,Yrotateplot3];
    
    Xplot4 = plot4(1,1);
    Yplot4 = plot4(1,2);
    angplot4 = rotateangle; 
    Xcplot4 = xrotatecenter;
    Ycplot4 = yrotatecenter;
    Xrotateplot4 =  (Xplot4-Xcplot4)*cosd(angplot4) + (Yplot4-Ycplot4)*sind(angplot4) + Xcplot4;
    Yrotateplot4 = -(Xplot4-Xcplot4)*sind(angplot4) + (Yplot4-Ycplot4)*cosd(angplot4) + Ycplot4;
    plot4rotate = [Xrotateplot4,Yrotateplot4];
    
    XplotC = Xcross;
    YplotC = Ycross;
    angplotC = rotateangle;
    XcplotC = xrotatecenter;
    YcplotC = yrotatecenter;
    XrotateplotC =  (XplotC-XcplotC)*cosd(angplotC) + (YplotC-YcplotC)*sind(angplotC) + XcplotC;
    YrotateplotC = -(XplotC-XcplotC)*sind(angplotC) + (YplotC-YcplotC)*cosd(angplotC) + YcplotC;
    plotCrotate = [XrotateplotC,YrotateplotC];

    rotateAroundd = ~rotateAroundd;
    filename1=[file_name_all_change,'a.jpg'];;
    imwrite(rotateAroundd,['F:\XXX\',filename1]);
    imwrite(BW,['F:\XXX\',filename1]);
 
    %% Reads the contour of anterior margin.
    
    maxdot_zhi = min(Yrotateplot3,Yrotateplot1);
    [row_zhi,column_zhi]=find(dotorigin_origin_rotate(:,1)<max(Xrotateplot3,Xrotateplot1)+300&dotorigin_origin_rotate(:,1)>min(Xrotateplot3,Xrotateplot1)-300&dotorigin_origin_rotate(:,2)<maxdot_zhi);
    dot_zhi = dotorigin_origin_rotate(row_zhi,:);
    dot_bottom = dotorigin_origin_rotate(find(dotorigin_origin_rotate(:,2)==max(dotorigin_origin_rotate(:,2))),:);
    dot_bottom(1,1)= XrotateplotC;
    dot_zhi(:,1) = dot_zhi(:,1)-dot_bottom(1,1);
    dot_zhi(:,2) = -(dot_zhi(:,2)-dot_bottom(1,2));
    dot_zhi = unique(dot_zhi,'row','stable');
    windowWidth = 45;
    polynomialOrder = 2;
    p_zhi = dot_zhi;
    name_zhi = [filename_zhi];               
    eval([[name_zhi,'= p_zhi'],';']); 
    %% plot
    subplot(1,2,1);
    scatter(p_zhi(:,1),p_zhi(:,2));
    save(['F:\XXX\', filename_zhi], [filename_zhi]);
   
    %% Reads the contour of heel region.
    mindot_gen = max(Yrotateplot4,Yrotateplot2);
    [row_gen,column_gen]=find(dotorigin_origin_rotate(:,1)<max(Xrotateplot4,Xrotateplot2)+300&dotorigin_origin_rotate(:,1)>min(Xrotateplot4,Xrotateplot2)-300&dotorigin_origin_rotate(:,2)>mindot_gen);
    dot_gen=dotorigin_origin_rotate(row_gen,:);
    dot_gen(:,1) = dot_gen(:,1)-dot_bottom(1,1);
    dot_gen(:,2) = -(dot_gen(:,2)-dot_bottom(1,2));
    dot_gen = unique(dot_gen,'row','stable');
 
    %% plot
    subplot(1,2,2);
    scatter(p_gen(:,1),-p_gen(:,2));
    filename2=[file_name_all_change,'b.jpg'];
    saveas(gca,['F:\XXX\',filename2]);
    save(['F:\XXX\', filename_gen], [filename_gen]);   
    clearvars -except List_all len_all iall
    
end