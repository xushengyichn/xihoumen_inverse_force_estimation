 function [S_a, S_v, S_d, n_sensors] = sensor_selection(loc_acc, loc_vel, loc_dis, node_loc, phi, nodeondeck, Mapping_data)

        if ~isempty(loc_acc)
            acc_node = FindNodewithLocation(loc_acc, node_loc, nodeondeck);
            acc_node_list = reshape(permute(acc_node, [2 1]), [], 1); % 交错重塑
            % acc_node_list=acc_node(:);
        else
            acc_node_list = [];
        end

        if ~isempty(loc_vel)
            vel_node = FindNodewithLocation(loc_vel, node_loc, nodeondeck);
            vel_node_list = reshape(permute(vel_node, [2 1]), [], 1);
        else
            vel_node_list = [];
        end

        if ~isempty(loc_dis)
            dis_node = FindNodewithLocation(loc_dis, node_loc, nodeondeck);
            dis_node_list = reshape(permute(dis_node, [2 1]), [], 1);
        else
            dis_node_list = [];
        end

        n_acc = length(acc_node_list);
        n_vel = length(vel_node_list);
        n_dis = length(dis_node_list);
        acc_matrix_seq = node2matrixseq(acc_node_list, Mapping_data);
        vel_matrix_seq = node2matrixseq(vel_node_list, Mapping_data);
        dis_matrix_seq = node2matrixseq(dis_node_list, Mapping_data);

        n_sensors = n_acc + n_vel + n_dis;

        S_a = zeros(n_sensors, size(phi, 1));
        S_v = zeros(n_sensors, size(phi, 1));
        S_d = zeros(n_sensors, size(phi, 1));

        for k1 = 1:n_sensors

            if k1 <= n_acc
                S_a(k1, acc_matrix_seq(k1)) = 1;
            elseif k1 <= n_acc + n_vel
                S_v(k1, vel_matrix_seq(k1 - n_acc)) = 1;
            else
                S_d(k1, dis_matrix_seq(k1 - n_acc - n_vel)) = 1;
            end

        end

    end
