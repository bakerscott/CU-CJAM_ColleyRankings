% APPM Matrix Methods Project
% Robert Hakulin
% Colley Ranking

clc; close all; clear all;

%% Sets up the C or Colley matrix 
% Loads in data 
data = load('4Data.csv');
data = data(1:2:end,:);
fid = fopen('4Teams.csv');
T = textscan(fid,'%s','delimiter',','); T = T{1};
fclose(fid);

% Extracts team names and corresponding team number from cell array
for i = 1:max(size(T))/2
    tnames{i} = T{2*i-1};
    tnums(i) = str2num(T{2*i});
end

% Initializes Colley Matrix and b vector
C = diag(2*ones(max(size(tnums)),1));
b = ones(max(size(tnums)),1);
wcnt = zeros(1,max(size(tnums)));
lcnt = zeros(1,max(size(tnums)));

% Loops through data and updates Colley and b vector
for i = 1:size(data)
    t1 = data(i,1); t2 = data(i,2);
    wol = data(i,3);
    
    % Updates off diagnols
    C(t1,t2) = C(t1,t2)-1;
    C(t2,t1) = C(t2,t1)-1;
    
    % Updates diagnol
    C(t1,t1) = C(t1,t1)+1;
    C(t2,t2) = C(t2,t2)+1;
    
    % Updates win/lose counter
    if wol == 1
        wcnt(t1) = wcnt(t1)+1;
        lcnt(t2) = lcnt(t2)+1;
    else
        wcnt(t2) = wcnt(t2)+1;
        lcnt(t1) = lcnt(t1)+1;
    end
end

% Loops through b vector
for j = 1:length(b)
    b(j) = b(j)+(wcnt(j)-lcnt(j))/2;
end

% Solves for rating vector, r from Cr=b
r = C\b;

% Sorts ranking team name, and team number vectors from largest to smallest
order = sort(r,'descend');
for k = 1:length(r)
    cnt = 1;
    while r(k) ~= order(k)
        if cnt ~= length(r)
            if r(cnt) <= r(cnt+1)
                temp = r(cnt);
                temp2{1} = tnames{cnt};
                r(cnt) = r(cnt+1);
                r(cnt+1) = temp;
                tnames{cnt} = tnames{cnt+1};
                tnames{cnt+1} = temp2{1};
            end
            cnt = cnt+1;
        else
            cnt = 1;
        end
    end
end

% Tabulates the ordered results
Rating = r;
Rank = (1:1:length(r))';
Results = table(Rank,Rating,'RowNames',tnames)
            