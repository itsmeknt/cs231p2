% K means clustering

function [assignment, centers]=kMeans(inputData,k)

    % kMeans accepts as input a set of m data vectors (as rows), each vector
    % has length n. It also accepts k, the number of clusters. It returns a 
    % vector, assignments, in which assignment(j) gives the cluster number 
    % to which data vector j is assigned. It also returns the coordinates 
    % of each cluster center in centers, i.e. center(i,:) gives the center of
    % cluster i.

    [m,n]=size(inputData);

    % Initialize cluster centers
    minData=min(inputData);
    maxData=max(inputData);
    interval=(maxData-minData)/k;
    centers=zeros(k,n);
    centers_new=centers;
    for i=1:k
        centers(i,:)=minData+interval/2;
        minData=minData+interval;
    end
    
    assignment=zeros(m,1);
    converged=0;
    
    % iterate until convergence
    while (converged~=1)
        % Assign each data point to the nearest cluster center
        for i=1:m
           assignment(i,1)=getClusterCenter(centers,inputData(i,:));
        end
        % Update cluster centers
        for j=1:k
            points=find(assignment==j);
            if (numel(points)~=0)
                centers_new(j,:)=mean(inputData(points,:),1);
            end
        end
        % Check for convergence
        if (norm(centers_new-centers)<1)
            converged=1;
            centers=centers_new;
        else
            centers=centers_new;
        end
    end
end

function assignmt=getClusterCenter(centers,datapoint)
% Assigns datapoint to its  cluster

k=size(centers,1);
distVectors=centers-repmat(datapoint,k,1);
distances=sum(distVectors.^2,2);
[~,idx]=min(distances);
assignmt=idx;
end
