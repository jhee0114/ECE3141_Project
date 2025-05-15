%% Simulating Hash Table Performance

% Making the hash function
function index = simpleHash(ip, tableSize)
    %prime number used for hashing
    prime = 31; 
    hash = mod(prime*double(ip), 2^32); %Using 2^32 as int is 32 bits to avoid overflow

    index = mod(hash, tableSize)+1; %matlab indexing starts at 1
end

% Insert function for hash table. With chaining in case of collision
function hashTable = insert(hashTable, ip, next, tableSize)
    index = simpleHash(ip, tableSize);
    
    %Check if the cell is empty
    if isempty(hashTable{index})
        hashTable{index} = {ip, next};
    else
        hashTable{index} = [hashTable{index}; {ip, next}];
    end
end

% Lookup function for the hash table
function result = lookup(hashTable, ip, tableSize)
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


%Generate a list of IP addresses
function ipList = generateRandomIP(tableSize)
    ipList = randi([0, 2^32-1], tableSize, 1, 'uint32');
end

tableSize = 100;
hashTable = cell(1, tableSize);

%Generate the hash table of ips
ipList = generateRandomIP(tableSize);
for i = 1:length(ipList)
    ip = ipList(i);
    insert(hashTable, ip, i, tableSize);
end

%Time lookup time
tic;
for i = 1:length(ipList)
    ip = ipList(i);
    next = lookup(hashTable, ip, tableSize);
end
elapsedTime = toc;
fprintf('Average lookup time: %.6f ms\n', (elapsedTime / tableSize) * 1000);



