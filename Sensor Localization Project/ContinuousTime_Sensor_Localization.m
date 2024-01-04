clear
clc
anchors = [0 0; 10 0; 5 10]; %setting anchor points

m = 2; %2D
node=10;

identity = eye(node);
initialValue = abs(10*rand(node, 1));
%initialValue = abs(10*rand(nodes, 1))

l = 0;
r = 16;
sidelength = r/2;
points=locate(anchors, m, node);

m1 = [anchors(:,2); anchors(1,2)];
m2 = [anchors(:,1); anchors(1,1)];

%intital value for each axis
x=initialValue;
y=initialValue;
initCond = initialConditions(node, 500);

%finding close node... this may need to change to get better accuracy
close_node = sum(distance_between(anchors, points, node)<(sidelength));
minimum = close_node(1:node);
while min(close_node(1:node)) <(m+2) %while min(minimum) <(m+2) 
    l=l+1;
   
    points = locate(anchors, m, node);
    
    close_node = sum(distance_between(anchors, points, node)<(sidelength));
end

[itr, A, B] = update_location(anchors, r, m, points, node);
        % = update_location(anchors, r, m, points, node);
        %A = update_location(anchors, r, m, points, node);
        %B = update_location(anchors, r, m, points, node);
    l=0;
    while itr==false
        l=l+1;
        points = locate(anchors, m, node);

        
        initCond = initialConditions(node, 500);
        
        close_node = sum(distance_between(anchors, points, node)<(sidelength));
        minimum = min(close_node(1:node));
        while min(close_node(1:node))<(m+2)
            l=l+1;
           
            points = locate(anchors, m, node);
            
            close_node = sum(distance_between(anchors, points, node)<(sidelength));
        end

        [itr, A, B] = update_location(anchors, r, m, points, node);
        %A = update_location(anchors, r, m, points, node);
        %B = update_location(anchors, r, m, points, node);
        
    end

%syms t;
t=0:0.01:10; %greater time step means smoother graph
%%sys_a =ss(A,B,identity,[]);
%for i=1:node
       
        %ct
        %l = 0
        %if lastpoint == nextpoint
            %l=1;
            
            %while(length(prediction)>m+1) 
                %lsim(sys_a)
                %
replicatey = repmat(anchors(:,2), 1, length(t));
replicatex = repmat(anchors(:,1), 1, length(t));

sys_a = ss(A-identity,B,identity,[]);
%sys_b =ss(A,B,identity,[]);
%initial conditions
%dx=lsim(sys_a, replicatex, t, x)';
%dy=lsim(sys_a, replicatey, t, y)';




%plotting
dx=(lsim(sys_a, replicatex, t, x))';
dy=(lsim(sys_a, replicatey, t, y))';


%actual values 


figure(); 
plot(m2,m1)
hold on


plot(points(:,1), points(:,2),'x','Color','r','Linewidth',1)
%for l = 1:len(node)
for l = 1:node
    
    %graphing sensors
    plot(dx(l,:),dy(l,:),'-')


    
end
%clear


function dist_to_node = distance_between(anchors, points, node)
    emptyArray = zeros(node+length(anchors));
    dist_to_node = emptyArray;
    for I=1:node+length(anchors)
        for l=1:node+length(anchors)
            if I>node && l>node
                dist_to_node(I,l) = sqrt(sum((anchors(I-node,:)-anchors(l-node,:)).^2));
            else
                if l>node
                    ctdistance = sqrt(sum((points(I,:)-anchors(l-node,:)).^2));
                    dist_to_node(I,l) = ctdistance;
                end

                if I>node
                    ctdistance = sqrt(sum((anchors(I-node,:)-points(l,:)).^2));
                    dist_to_node(I,l) = ctdistance; 
                end
                
            if I<node && l<node
                dist_to_node(I,l)=sqrt(sum((points(I,:)-points(l,:)).^2));
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
    initialarray=zeros(1,m+1);
    for i=1:node
        initialarray(1) = rand(1);
        for l=2:m+1

            initialarray(l)=(1-sum(initialarray(1:(l-1))))*rand(1); %fill with random num
        end


        initialarray=initialarray(randperm(m+1));
        for j=1:m
            points(i,j)=initialarray * anchors(:,j);
        end
    end
end


function [A,B]=get_old_location(anchors, prediction, points, m, i)
    designated_pts = [points; anchors];
    
    %default setting
    A=zeros(1,length(points));
    B=zeros(1,m+1);
    size = length(prediction);
    amountofpts = length(points);
    
    for l = 1 : size
       %adjust the localization of the pts
       node=m==sum(designated_pts == prediction(l, :),2);
       nodeNew=find(node);
       adjustedpred = [prediction([1:l-1,l+1:end],:);points(i,:)];
       [~,sizeof] = convhull(prediction);
       
       [~,actualsizeof] = convhull(adjustedpred);
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
    sidelength = r/2;
    itr=false;

    A=zeros(node);
    B=zeros(node,m+1);

 
    distance=distance_between(anchors, points, node)<(sidelength);


    for i=1:node
       
        distancearray=distance(:,i);
        distancearray(i)=0;
        ctprediction1 = points(distancearray(1:node),:); 
        ctprediction2 = anchors(distancearray(1+node:end),:);
        prediction=[ctprediction1; ctprediction2];



        newmat = [prediction;points(i,:)];
        [~, convex_hull]=convhull(prediction);
        [~, convex_hull_mod]=convhull(newmat);
        
        if convex_hull==convex_hull_mod
            l=1;
            
            while(length(prediction)>m+1) 
                current=prediction([1:l-1,1+l:end],:);
                currentpointgroup = [current; points(i,:)];
                [~, convex_hull_mod]=convhull(currentpointgroup);


                
                [~, convex_hull] = convhull(current);
                
                if convex_hull == convex_hull_mod %check if they have same value
                    prediction = current; %location guessed
                    l = 1;
                else
                    l = l + 1;
                end
            end
        else 
            itr=false;
            return 
        end
            [A(i,:), B(i,:)] = get_old_location(anchors, prediction, points, m, i);
            
            

    end
    itr=true;
end


