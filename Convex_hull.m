[k,av] = convhull(dotorigin_origin);
if length(dotorigin_origin_bottom)~=length(dotorigin_origin)
    [k_bottom,av_bottom] = convhull(dotorigin_origin_bottom);
    dot_add=[dotorigin_origin(k,:);dotorigin_origin_bottom(k_bottom,:)];
    [k,av] = convhull(dot_add);
    dotconvhull = [dot_add(k,1),dot_add(k,2)];
    s = length(dotconvhull);
    for ii=1:s-1
        convhulldist(ii,1) = sqrt((dotconvhull(ii,1)-dotconvhull(ii+1,1))^2+(dotconvhull(ii,2)-dotconvhull(ii+1,2))^2);
    end
    maxvalue1 = max(convhulldist);
    max1 = find(convhulldist==maxvalue1);
    convhulldist(max1,1) = 0;
    maxvalue2 = max(convhulldist);
    max2 = find(convhulldist==maxvalue2);
    temppoint(1,:) = dot_add(k(max1,1),:);
    temppoint(2,:) = dot_add(k(max1+1,1),:);
    temppoint(3,:) = dot_add(k(max2+1,1),:);
    temppoint(4,:) = dot_add(k(max2,1),:);
    tempdist(1,1) = abs(temppoint(1,2)-temppoint(2,2));
    tempdist(2,1) = abs(temppoint(1,2)-temppoint(3,2));
    tempdist(3,1) = abs(temppoint(1,2)-temppoint(4,2));
    mintempdist = min(tempdist);
    min1 = find(tempdist==mintempdist)+1;
    if temppoint(1,1)-temppoint(min1,1)<0
        plot1 = temppoint(1,:);
        plot3 = temppoint(min1,:);
    else
        plot3 = temppoint(1,:);
        plot1 = temppoint(min1,:);
    end
    temppoint_deleete = temppoint;
    temppoint_deleete(min1,:)=[];
    temppoint_deleete(1,:)=[];
    if temppoint_deleete(1,1)-temppoint_deleete(2,1)<0
        plot2 = temppoint_deleete(1,:);
        plot4 = temppoint_deleete(2,:);
    else
        plot4 = temppoint_deleete(1,:);
        plot2 = temppoint_deleete(2,:);
    end
    if abs(plot1(1,1)-plot3(1,1))>abs(plot2(1,1)-plot4(1,1))
    else
        temp1=plot1;
        temp3=plot3;
        plot1=plot2;
        plot3=plot4;
        plot2=temp1;
        plot4=temp3;
    end
    dotorigin_origin=[dotorigin_origin;dotorigin_origin_bottom];
else
    dotconvhull = [dotorigin_origin(k,1),dotorigin_origin(k,2)];
    s = length(dotconvhull);
    for ii=1:s-1
        convhulldist(ii,1) = sqrt((dotconvhull(ii,1)-dotconvhull(ii+1,1))^2+(dotconvhull(ii,2)-dotconvhull(ii+1,2))^2);
    end
    maxvalue1 = max(convhulldist);
    max1 = find(convhulldist==maxvalue1);
    convhulldist(max1,1) = 0;
    maxvalue2 = max(convhulldist);
    max2 = find(convhulldist==maxvalue2);
    temppoint(1,:) = dotorigin_origin(k(max1,1),:);
    temppoint(2,:) = dotorigin_origin(k(max1+1,1),:);
    temppoint(3,:) = dotorigin_origin(k(max2+1,1),:);
    temppoint(4,:) = dotorigin_origin(k(max2,1),:);
    tempdist(1,1) = abs(temppoint(1,2)-temppoint(2,2));
    tempdist(2,1) = abs(temppoint(1,2)-temppoint(3,2));
    tempdist(3,1) = abs(temppoint(1,2)-temppoint(4,2));
    mintempdist = min(tempdist);
    min1 = find(tempdist==mintempdist)+1;
    if temppoint(1,1)-temppoint(min1,1)<0
        plot1 = temppoint(1,:);
        plot3 = temppoint(min1,:);
    else
        plot3 = temppoint(1,:);
        plot1 = temppoint(min1,:);
    end
    temppoint_deleete = temppoint;
    temppoint_deleete(min1,:)=[];
    temppoint_deleete(1,:)=[];
    if temppoint_deleete(1,1)-temppoint_deleete(2,1)<0
        plot2 = temppoint_deleete(1,:);
        plot4 = temppoint_deleete(2,:);
    else
        plot4 = temppoint_deleete(1,:);
        plot2 = temppoint_deleete(2,:);
    end
    if abs(plot1(1,1)-plot3(1,1))>abs(plot2(1,1)-plot4(1,1))
    else
        temp1=plot1;
        temp3=plot3;
        plot1=plot2;
        plot3=plot4;
        plot2=temp1;
        plot4=temp3;
    end
end
