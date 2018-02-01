% APPM Matrix Methods Project
% Robert Hakulin
% ELO Ranking

clc; close all; clear all;

%% Sets up the C or Colley matrix 
% Loads in data 
data = load('Div1Data.csv');
%data = data(1:2:end,:);
fid = fopen('Div1Teams.csv');
T = textscan(fid,'%s','delimiter',','); T = T{1};
fclose(fid);

% Extracts team names and corresponding team number from cell array
for i = 1:max(size(T))/2
    tnames{i} = T{2*i-1};
    tnums(i) = str2num(T{2*i});
end

% Initializes Colley Matrix and b vector
C = diag(2*ones(max(size(tnums)),1));
b = 1000*ones(max(size(tnums)),1);
k = 50;

% Loops through data and updates Colley and b vector
for i = 1:size(data)
    t1 = data(i,1); t2 = data(i,2);
    wol = data(i,3);
    p1c = 10^(b(t1)/400)/(10^(b(t2)/400)+10^(b(t1)/400));
    p2c = 10^(b(t2)/400)/(10^(b(t2)/400)+10^(b(t1)/400));
    b(t1) = b(t1)+k*(wol-p1c);
    b(t2) = b(t2)+k*(abs(wol-1)-p2c);
end


% Solves for rating vector, r from Cr=b
r = b;

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
ELO = round(r);
Rank = (1:1:length(r))';
Results = table(Rank,ELO,'RowNames',tnames)
            