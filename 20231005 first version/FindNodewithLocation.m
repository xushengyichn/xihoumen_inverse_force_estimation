    function node = FindNodewithLocation(loc, node_loc, nodeondeck)
        %myFun - Description
        %
        % Syntax: result = findnode(input)
        %
        % Long description
        node_seq = zeros(1, length(loc));

        for k1 = 1:length(loc)
            node_seq(k1) = find(node_loc >= loc(k1), 1);
        end

        node = nodeondeck(node_seq, :);
    end
