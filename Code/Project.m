clear all; close all; clc;
%% Normal Hash table

% Making the hash function
function index = simpleHash(ip, tableSize)
    %prime number used for hashing
    prime = 31; 
    hash = mod(prime*double(ip), 2^32); %Using 2^32 as int is 32 bits to avoid overflow

    index = mod(hash, tableSize)+1; %matlab indexing starts at 1
end

% Insert function for hash table. With chaining in case of collision
function [hashTable, collisionCount] = hashInsert(hashTable, ip, next, tableSize, collisionCount)
    index = simpleHash(ip, tableSize);
    
    %Check if the cell is empty
    if isempty(hashTable{index})
        hashTable{index} = {ip, next};
    else
        hashTable{index} = [hashTable{index}; {ip, next}];
        collisionCount = collisionCount + 1;
    end
end

% Lookup function for the hash table
function result = hashLookup(hashTable, ip, tableSize)
    index = simpleHash(ip, tableSize);
    bucket = hashTable{index}; %bucket because of chaining

    result = "Miss"; %Assume it is a miss at the beginning
    if ~isempty(bucket)
        for i = 1:size(bucket, 1)
            if isequal(bucket{i, 1}, ip)
                result = bucket{i, 2};
                return;
            end
        end
    end
end


%% D-Left Hashing

% Hash function with a seed
function index = dLeftHash(ip, subtableSize, seed)
    prime = 31 + seed*2;
    hash = mod(prime*double(ip), 2^32);
    index = mod(hash, subtableSize)+1;
end

% Insert function for d-left hashing
function [hashTable, collisionCount] = dLeftInsert(hashTable, ip, next, d, subtableSize, collisionCount)
    minLen = inf; %Starting off at infinity for comparison
    chosenSub = 1; %Default starting sub table
    
    %Checking the collision
    isCollision = true;
    for i = 1:d
        index = dLeftHash(ip, subtableSize, i);
        if isempty(hashTable{i}{index}) %check if that space is empty
            hashTable{i}{index} = struct('ip', ip, 'value', next);
            isCollision = false;
            return;
        else 
            currentLen = length(hashTable{i}{index}); %Check the current length to find the subtable with the shortest length
            if currentLen < minLen
                minLen = currentLen;
                chosenSub = i;
                bestIndex = index;
            end
        end
    end

    %Insert the IP address into the least loaded subTable
    hashTable{chosenSub}{bestIndex}(end+1) = struct('ip', ip, 'value', next);

    %Check if collision has occurred or not
    if isCollision
        collisionCount = collisionCount + 1;
    end
end


%Look up function for d-left hashing
function result = dLeftLookup(hashTable, ip, d, subtableSize)
    result = "Miss";
    for i = 1:d %Look through all of the sub tables
        index = dLeftHash(ip, subtableSize, i);
        bucket = hashTable{i}{index};
        if ~isempty(bucket) % If the bucket is not empty
            for j = 1:length(bucket) %Look through all of the chained values
                if bucket(j).ip==ip
                    result = bucket(j).value;
                    return;
                end
            end
        end
    end
    %Otherwise, it is a miss
end



%% Trie-Based Lookup

%Function to make a trie node
function node = createTrieNode()
    node.left = []; %For bit 0
    node.right = []; %For bit 1
    node.value = []; %The ip value
end

%Convert the ip to a binary string to create the tree
function binStr = ipToBinary(ip)
    binStr = dec2bin(ip, 32); %Return the 32 bit binary string
end

% Insert the ip into the tree
function trie = trieInsert(trie, ip, value)
    binStr = ipToBinary(ip);
    node = trie;
    
    % Loop through all of the characters in the ip address
    for i = 1:32
        % If the binary is 0, append to the left
        if binStr(i) == "0"
            if isempty(node.left)
                node.left = createTrieNode();
            end
            node = node.left;
        % IF the binary is 1, append to the right
        else
            if isempty(node.right)
                node.right = createTrieNode();
            end
            node = node.right;
        end
    end
    node.value = value;
end

%Traverse through the tree
%function result = trieLookup(trie, ip)
%    binStr = ipToBinary(ip);
%    node = trie;
%    result = "Miss";
    
%    %Traverse through all of the values in the binary string
%    for i = 1:32
%        if isempty(node)
%            return;
%        end
%        if ~isempty(node.value)
%            result = node.value;
%        end
%        if (binStr(i))=='0'
%            node = node.left;
%        else
%            node = node.right;
%        end
%    end
    
%    %Return the value of the node
%    if ~isempty(node) && ~isempty(node.value)
%        result = node.value;
%    end
%end
function result = trieLookup(trie, ip)
    binStr = ipToBinary(ip);
    node = trie;

    for i = 1:32
        if isempty(node)
            result = "Miss";
            return;
        end
        if binStr(i) == '0'
            node = node.left;
        else
            node = node.right;
        end
    end

    % After finishing all 32 bits, check if exact match is found
    if ~isempty(node) && ~isempty(node.value)
        result = node.value;
    else
        result = "Miss";
    end
end


%% Testing All the different algorithms

%Generate a list of IP addresses
function ipList = generateRandomIP(numberOfIPs)
    ipList = randi([0, 2^32-1], numberOfIPs, 1, 'uint32');
end

ipCounts = 1000:1000:10000;
hashLookupTimes = zeros(size(ipCounts));
dleftLookupTimes = zeros(size(ipCounts));
trieLookupTimes = zeros(size(ipCounts));

hashCollisions = zeros(size(ipCounts));
dleftCollisions = zeros(size(ipCounts));

hashInsertTimes = zeros(size(ipCounts));
dleftInsertTimes = zeros(size(ipCounts));
trieInsertTimes = zeros(size(ipCounts));

for idx = 1:length(ipCounts)
    numberOfIPs = ipCounts(idx);

    ipList = generateRandomIP(numberOfIPs);
    


    %--- Hash Table----
    %Generate the hash table of ips
    collisionCountHash = 0;
    hashTable = cell(1, numberOfIPs);
    tic;
    for i = 1:length(ipList)
        ip = ipList(i);
        [hashTable, collisionCountHash] = hashInsert(hashTable, ip, i, numberOfIPs, collisionCountHash);
    end
    elapsedTimeHash = toc;
    hashInsertTimes(idx) = elapsedTimeHash/numberOfIPs;
    hashCollisions(idx) = collisionCountHash;

    %Time lookup time for Hash Table
    tic;
    for i = 1:length(ipList)
        ip = ipList(i);
        next = hashLookup(hashTable, ip, numberOfIPs);
    end
    elapsedTimeHash = toc;
    hashLookupTimes(idx) = elapsedTimeHash/numberOfIPs;



     %--- d-left Hash Table ----
    
    %Generate the d-left hash table of ips
    d = 4; %number of sub tables
    subtableSize = ceil(numberOfIPs/4); %Must contain all of the entries
    
    %Iniitilize d-left hash table
    collisionCountdLeft = 0;
    dLeftHashTable = cell(1, d);
    for i = 1:d
        dLeftHashTable{i} = cell(1, subtableSize);
    end
    
    %Insert the IPs into the d-left hashTable
    tic;
    for i = 1:length(ipList)
        ip = ipList(i);
        [dLeftHashTable, collisionCountdLeft] = dLeftInsert(dLeftHashTable, ip, i, d, subtableSize, collisionCountdLeft);
    end
    elapsedTimedLeft = toc;
    dleftCollisions(idx) = collisionCountdLeft;
    dleftInsertTimes(idx) = elapsedTimedLeft/numberOfIPs;

    %Time lookup time for d-Left Hash Table
    tic;
    for i = 1:length(ipList)
        ip = ipList(i);
        next = dLeftLookup(dLeftHashTable, ip, d, subtableSize);
    end
    elapsedTimedLeft = toc;
    dleftLookupTimes(idx) = elapsedTimedLeft/numberOfIPs;



    % --- Trie ---
    trie = createTrieNode();
    
    %Generate the Trie 
    tic
    for i = 1:numberOfIPs
        trie = trieInsert(trie, ipList(i), i);
    end
    elapsedTimeTrie = toc;
    trieInsertTimes(idx) = elapsedTimeTrie/numberOfIPs;
    
    % Time the lookup
    tic;
    for i = 1:numberOfIPs
        result = trieLookup(trie, ipList(i));
    end
    elapsedTimeTrie = toc;
    trieLookupTimes(idx) = elapsedTimeTrie/numberOfIPs;
end

% Plotting the average lookup time results
figure(1);
plot(ipCounts, hashLookupTimes * 1000, '-', 'LineWidth', 2);
hold on;
plot(ipCounts, dleftLookupTimes * 1000, '-', 'LineWidth', 2);
plot(ipCounts, trieLookupTimes * 1000, '-', 'LineWidth', 2);
hold off;
xlabel('Number of IPs');
ylabel('Average Lookup Time (ms)');
title('Average Lookup Times: Hash vs D-Left vs Trie');
legend('Hash Table', 'D-Left Hashing', 'Trie','Location','best');
grid on;

%Plotting the average insert time results
figure(2);
plot(ipCounts, hashInsertTimes * 1000, '-', 'LineWidth', 2);
hold on;
plot(ipCounts, dleftInsertTimes * 1000, '-', 'LineWidth', 2);
plot(ipCounts, trieInsertTimes * 1000, '-', 'LineWidth', 2);
hold off;
xlabel('Number of IPs');
ylabel('Average Insert Time (ms)');
title('Average Insert Times: Hash vs D-Left vs Trie');
legend('Hash Table', 'D-Left Hashing', 'Trie','Location','best');
grid on;

%Plotting the number of collisions results
figure(3);
plot(ipCounts, hashCollisions, '-', 'LineWidth', 2);
hold on;
plot(ipCounts, dleftCollisions, '-', 'LineWidth', 2);
hold off;
xlabel('Number of IPs');
ylabel('Number of Collisions)');
title('Number of Collisions: Hash vs D-Left');
legend('Hash Table', 'D-Left Hashing','Location','best');
grid on;
