   
    function matrix_seq = node2matrixseq(node_list, KMmapping)
        UYNode = KMmapping.Node(KMmapping.DOF == 'UY');
        UYMatrix = KMmapping.MatrixEqn(KMmapping.DOF == 'UY');
        indices = [];

        for k1 = 1:length(node_list)
            indices(k1) = find(UYNode == node_list(k1));
        end

        if isempty(indices)
            matrix_seq = [];
        else
            matrix_seq = UYMatrix(indices);
        end

    end