%setting anchor points

anchors = [0 0; 10 0; 5 10];
m = 2;
l=0;
iterations=1000;





%initialValue = abs(10*rand(nodes, 1))
r = 15;


%to choose how many nodes you want, change the value of 'node'
%FOR m = 4
node = 4;
sidelength = r/2



%itr = true; %may need to remove this

points = locate(anchors, m, node);


initCond = initialConditions(node, 500)

close_node=sum(distance_between(anchors, points, node)<(sidelength));
miniumdist = min(close_node(1:node));


while miniumdist<(m+2)
    l=l+1;
   
    points = locate(anchors, m, node);
    
    close_node = sum(distance_between(anchors, points, node)<(sidelength));
end


[itr, A, B]=update_location(anchors, r, m, points, node);





l=0;
while itr==false
    l=l+1;
   points = locate(anchors, m, node);
    l=0;
   
    close_node = sum(distance_between(anchors, points, node)<(sidelength));
    miniumdist = min(close_node(1:node));
    while miniumdist<(m+2)
        l=l+1;
       
        points = locate(anchors, m, node);
        

    end

    [itr, A, B] = update_location(anchors, r, m, points, node);
   
end




xaxis=initCond;
yaxis=initCond;

for l=2:iterations
    xaxis(:,l)=A*xaxis(:,l-1)+B*anchors(:,1);
    yaxis(:,l)=A*yaxis(:,l-1)+B*anchors(:,2);
end

figure();

plot([anchors(:,1);anchors(1,1)],[anchors(:,2);anchors(1,2)],'-x')
hold on

plot(points(:,1),points(:,2),'x','Color','r','Linewidth',3)

for l = 1:node
    plot(xaxis(l,:), yaxis(l,:)) %might be beter than


end
clc

%%%%FOR M=10

%setting anchor points

anchors = [0 0; 10 0; 5 10];
m = 2;
l=0;
iterations=1000;





%initialValue = abs(10*rand(nodes, 1))
r = 15;


%to choose how many nodes you want, change the value of 'node'
%FOR m = 10
node = 10;
sidelength = r/2



%itr = true; %may need to remove this

points = locate(anchors, m, node);


initCond = initialConditions(node, 500)

close_node=sum(distance_between(anchors, points, node)<(sidelength));
miniumdist = min(close_node(1:node));


while miniumdist<(m+2)
    l=l+1;
   
    points = locate(anchors, m, node);
    
    close_node = sum(distance_between(anchors, points, node)<(sidelength));
end


[itr, A, B]=update_location(anchors, r, m, points, node);





l=0;
while itr==false
    l=l+1;
   points = locate(anchors, m, node);
    l=0;
   
    close_node = sum(distance_between(anchors, points, node)<(sidelength));
    miniumdist = min(close_node(1:node));
    while miniumdist<(m+2)
        l=l+1;
       
        points = locate(anchors, m, node);
        

    end

    [itr, A, B] = update_location(anchors, r, m, points, node);
   
end




xaxis=initCond;
yaxis=initCond;

for l=2:iterations
    xaxis(:,l)=A*xaxis(:,l-1)+B*anchors(:,1);
    yaxis(:,l)=A*yaxis(:,l-1)+B*anchors(:,2);
end

figure();

plot([anchors(:,1);anchors(1,1)],[anchors(:,2);anchors(1,2)],'-x')
hold on

plot(points(:,1),points(:,2),'x','Color','r','Linewidth',3)

for l = 1:node
    plot(xaxis(l,:), yaxis(l,:)) %might be beter than


end

clc




%%%%FOR M=50

%setting anchor points

anchors = [0 0; 10 0; 5 10];
m = 2;
l=0;
iterations=1000;





%initialValue = abs(10*rand(nodes, 1))
r = 15;


%to choose how many nodes you want, change the value of 'node'
%FOR m = 10
node = 50;
sidelength = r/2



%itr = true; %may need to remove this

points = locate(anchors, m, node);


initCond = initialConditions(node, 500)

close_node=sum(distance_between(anchors, points, node)<(sidelength));
miniumdist = min(close_node(1:node));


while miniumdist<(m+2)
    l=l+1;
   
    points = locate(anchors, m, node);
    
    close_node = sum(distance_between(anchors, points, node)<(sidelength));
end


[itr, A, B]=update_location(anchors, r, m, points, node);





l=0;
while itr==false
    l=l+1;
   points = locate(anchors, m, node);
    l=0;
   
    close_node = sum(distance_between(anchors, points, node)<(sidelength));
    miniumdist = min(close_node(1:node));
    while miniumdist<(m+2)
        l=l+1;
       
        points = locate(anchors, m, node);
        

    end

    [itr, A, B] = update_location(anchors, r, m, points, node);
   
end




xaxis=initCond;
yaxis=initCond;

for l=2:iterations
    xaxis(:,l)=A*xaxis(:,l-1)+B*anchors(:,1);
    yaxis(:,l)=A*yaxis(:,l-1)+B*anchors(:,2);
end

figure();

plot([anchors(:,1);anchors(1,1)],[anchors(:,2);anchors(1,2)],'-x')
hold on

plot(points(:,1),points(:,2),'x','Color','r','Linewidth',3)

for l = 1:node
    plot(xaxis(l,:), yaxis(l,:)) %might be beter than


end

clc









%%%%%%%%%%%%%%%%%%%
%%%%%For R3
l=0;
m_r3=3;
nodes_r3=10;
rr3 = 15;
step=1000;
anchors_r3=[0 0 0; 10 0 0; 5 10 0; 5 5 10];

%generating the points

points = locate(anchors_r3, m_r3, nodes_r3);
initCond = initialConditions(nodes_r3, 500)


close_node=sum(distance_between(anchors_r3, points, nodes_r3)<(sidelength));
miniumdist = min(close_node(1:nodes_r3));


while min(close_node(1:nodes_r3))<(m_r3+2)
    l=l+1;
   
    points = locate(anchors_r3, m_r3, nodes_r3);
    
    close_node=sum(distance_between(anchors_r3, points, nodes_r3)<(sidelength));
end








[itr, A, B]=update_location(anchors_r3, rr3, m_r3, points, nodes_r3);

while itr==false
    l = l + 1;
    points = locate(anchors_r3, m_r3, nodes_r3);
    
    
    close_node=sum(distance_between(anchors_r3, points, nodes_r3)<(sidelength));
    miniumdist = min(close_node(1:nodes_r3));
    
    
    while min(close_node(1:nodes_r3))<(m_r3+2)
        l=l+1;
       
        points=locate(anchors_r3, m_r3, nodes_r3);
        
        close_node=sum(distance_between(anchors_r3, points, nodes_r3)<(sidelength));
    end
    [itr, A, B]=update_location(anchors_r3, rr3, m_r3, points, nodes_r3);

end



%check the spectral radius of A to ensure it is less than 1


%Initial conditions of x y and z
initCond = initialConditions(nodes_r3,50)
xaxis=initCond;
yaxis=initCond;
z=initCond;

%run the simulation
for l=5:step
    xaxis(:,l)=A*xaxis(:,l-1)+B*anchors_r3(:,1);
    yaxis(:,l)=A*yaxis(:,l-1)+B*anchors_r3(:,2);
    z(:,l)=A*z(:,l-1)+B*anchors_r3(:,3);
end

figure();

plot3([anchors_r3(:,1);anchors_r3(1,1)],[anchors_r3(:,2);anchors_r3(1,2)],[anchors_r3(:,3);anchors_r3(1,3)],'x')
hold on

plot3(points(:,1),points(:,2),points(:,3),'x','Color','r','Linewidth',1)

for l = 1:nodes_r3
    plot3(xaxis(l,:),yaxis(l,:),z(l,:),'-')
end





function dist_to_node = distance_between(anchors, points, node)
    emptyArray = zeros(node+length(anchors))
    dist_to_node = emptyArray
    for i=1:node+length(anchors)
        for j=1:node+length(anchors)
            if i>node && j>node
                dist_to_node(i,j)=sqrt(sum((anchors(i-node,:)-anchors(j-node,:)).^2));
            else
                if j>node
                    dist_to_node(i,j)=sqrt(sum((points(i,:)-anchors(j-node,:)).^2));    
                end

                if i>node
                    dist_to_node(i,j)=sqrt(sum((anchors(i-node,:)-points(j,:)).^2)); 
                end
                
            if i<node && j<node
                dist_to_node(i,j)=sqrt(sum((points(i,:)-points(j,:)).^2));
            end
            end
        end
    end
end

function initCond = initialConditions(nodes_r3, x)
    initCond=abs(10*rand(nodes_r3, x));
end

function points = locate(anchors, m, node)
    points=zeros(node,m);
    emptyarray=zeros(1,m+1);
    for i=1:node
        emptyarray(1) = rand(1);
        for l=2:m+1
            emptyarray(l)=(1-sum(emptyarray(1:(l-1))))*rand(1);
        end
        emptyarray=emptyarray(randperm(m+1));
        for j=1:m
            points(i,j)=emptyarray*anchors(:,j);
        end
    end
end


function [A,B]=get_old_location(anchors, prediction, points, m, i)
    designated_pts = [points; anchors];
    
    %default setting
    A=zeros(1,length(points));
    B=zeros(1,m+1);
    size = length(prediction)
    amountofpts = length(points)
    
    for l = 1 : size
       %adjust the localization of the pts
       node=m==sum(designated_pts == prediction(l, :),2);
       nodeNew=find(node);
       adjustedpred = [prediction([1:l-1,l+1:end],:);points(i,:)]
       [~,sizeof]=convhull(prediction);
       
       [~,actualsizeof]=convhull(adjustedpred);
       nodesize = actualsizeof / sizeof;

       if nodeNew<length(points)
           A(nodeNew) = nodesize;
       end
       
       if nodeNew>length(points)
           B(nodeNew-amountofpts) = nodesize;
       end
    end
end


function [itr, A, B] = update_location(anchors, r, m, points, node)
    sidelength = r/2
    itr=false;

    A=zeros(node);
    B=zeros(node,m+1);

    dist=distance_between(anchors, points, node);
    distance=dist<(sidelength);


    for i=1:node
       
        distancearray=distance(:,i);
        distancearray(i)=0;
        prediction=[points(distancearray(1:node),:);anchors(distancearray(node+1:end),:)];
        newmat = [prediction;points(i,:)];
        [~,convex_hull]=convhull(prediction);
        [~,convex_hull_mod]=convhull(newmat);
        
        if convex_hull==convex_hull_mod
            l=1;
            
            while(length(prediction)>m+1) 
                current=prediction([1:l-1,l+1:end],:);
                currentpointgroup = [current; points(i,:)];
                [~,convex_hull_mod]=convhull(currentpointgroup);



                [~,convex_hull] = convhull(current);
                
                if convex_hull == convex_hull_mod %check if they have same value
                    prediction = current; %location guessed
                    l = 1;
                else
                    l = l+1;
                end
            end
        else 
            itr=false;
            return 
        end
            [A(i,:),B(i,:)]=get_old_location(anchors, prediction, points, m, i);
            
            

    end
    itr=true;
end

