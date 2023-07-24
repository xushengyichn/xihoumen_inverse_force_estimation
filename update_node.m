function [node_modal_shape,nodes_displacement,nodes_displacementMagnitude,X_modal_shape,Y_modal_shape,Z_modal_shape]= update_node(node_modal_shape,modal_shape_plot,Mapping_data,elements)
    node_modal_shape_original = node_modal_shape;
    for k1= 1:length(modal_shape_plot)
        Node_temp= Mapping_data{k1,2};
        Node_DOF= Mapping_data{k1,3};
        if Node_DOF == "UX"
            Seq = find(node_modal_shape(:,1)==Node_temp);
            node_modal_shape(Seq,2) = node_modal_shape(Seq,2)+modal_shape_plot(k1);
        end
        if Node_DOF == "UY"
            Seq = find(node_modal_shape(:,1)==Node_temp);
            node_modal_shape(Seq,4) = node_modal_shape(Seq,4)+modal_shape_plot(k1);
        end%为了在matlab中绘图方便，这里注意一下，matlab中的Y轴对应的是Z轴，Z轴对应的是Y轴
        if Node_DOF == "UZ"
            Seq = find(node_modal_shape(:,1)==Node_temp);
            node_modal_shape(Seq,3) = node_modal_shape(Seq,3)+modal_shape_plot(k1);
        end
    
    end
    nodes_displacement = node_modal_shape(:,2:4)-node_modal_shape_original(:,2:4);
    nodes_displacementMagnitude = sqrt(sum(nodes_displacement.^2, 2));


    % Number of beams
    numBeams = size(elements, 1);

    for k1 = 1:numBeams
        % Get node indices for the i-th beam
        node1Index = elements{k1, 7};
        node2Index = elements{k1, 8};

        % Get the node coordinates
        node1(k1, :) = node_modal_shape(node_modal_shape(:, 1) == node1Index, 2:4);
        node2(k1, :) = node_modal_shape(node_modal_shape(:, 1) == node2Index, 2:4);
    end
    % Prepare matrix inputs for line function
    X_modal_shape = [node1(:,1), node2(:,1)].'; % X coordinates of line ends
    Y_modal_shape = [node1(:,2), node2(:,2)].'; % Y coordinates of line ends
    Z_modal_shape = [node1(:,3), node2(:,3)].'; % Z coordinates of line ends

end    