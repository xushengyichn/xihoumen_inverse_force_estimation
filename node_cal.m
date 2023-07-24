clc
clear
close all


node_loc = zeros(269,1);

for i = 2:11
   node_loc(i) = node_loc(i-1) + 5;
end

for i = 12:67
   node_loc(i) = node_loc(i-1) + 9;
end

for i = 68:83
   node_loc(i) = node_loc(i-1) + 3;
end

for i = 84:261
   node_loc(i) = node_loc(i-1) + 9;
end

for i = 262:269
   node_loc(i) = node_loc(i-1) + 3;
end

node_gap = diff(node_loc);
