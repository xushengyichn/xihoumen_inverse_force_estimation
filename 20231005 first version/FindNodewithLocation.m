    % function node = FindNodewithLocation(loc, node_loc, nodeondeck)
    %     %myFun - Description
    %     %
    %     % Syntax: result = findnode(input)
    %     %
    %     % Long description
    %     node_seq = zeros(1, length(loc));
    % 
    %     for k1 = 1:length(loc)
    %         node_seq(k1) = find(node_loc >= loc(k1), 1);
    %     end
    % 
    %     node = nodeondeck(node_seq, :);
    % end


    function nodes = FindNodewithLocation(loc, node_loc, nodeondeck)
    %myFun - Description
    %
    % Syntax: result = findnode(input)
    %
    % Long description
    nodes_seq_left = zeros(1, length(loc));
    nodes_seq_right = zeros(1, length(loc));
    
    for k1 = 1:length(loc)
        nodes_seq_right(k1) = find(node_loc >= loc(k1), 1);
        nodes_seq_left(k1) = max(1, nodes_seq_right(k1) - 1);

        distance_right = abs(node_loc(nodes_seq_right(k1)) - loc(k1));
        distance_left = abs(node_loc(nodes_seq_left(k1)) - loc(k1));
        total_distance = distance_right + distance_left;

        weights(1, k1) = distance_right / total_distance; % left node weight
        weights(2, k1) = distance_left / total_distance;  % right node weight
    end
    
    nodes.left = nodeondeck(nodes_seq_left, :);
    nodes.right = nodeondeck(nodes_seq_right, :);
    nodes.weights = weights;
end
