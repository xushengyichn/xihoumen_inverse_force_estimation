%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-10-12 22:19:31
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-11-01 20:29:10
%FilePath: \Exercises-for-Techniques-for-estimation-in-dynamics-systemsf:\git\xihoumen_inverse_force_estimation\20231005 first version\sensor_selection.m
%Description: 
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 function [S_a, S_v, S_d, n_sensors] = sensor_selection(loc_acc, loc_vel, loc_dis, node_loc, phi, nodeondeck, Mapping_data)

        if ~isempty(loc_acc)
            acc_node = FindNodewithLocation(loc_acc, node_loc, nodeondeck);
            acc_node_right = acc_node.right;
            acc_node_left = acc_node.left;
            acc_node_weights = acc_node.weights;
            acc_node_right_list = reshape(permute(acc_node_right, [2 1]), [], 1); % 交错重塑
            acc_node_left_list = reshape(permute(acc_node_left, [2 1]), [], 1); % 交错重塑
            % acc_node_list = reshape(permute(acc_node, [2 1]), [], 1); % 交错重塑
            % acc_node_list=acc_node(:);
        else
            % acc_node_list = [];
            acc_node_right_list = [];
            acc_node_left_list = [];
        end

        if ~isempty(loc_vel)
            vel_node = FindNodewithLocation(loc_vel, node_loc, nodeondeck);
            vel_node_right = vel_node.right;
            vel_node_left = vel_node.left;
            vel_node_weights = vel_node.weights;
            vel_node_right_list = reshape(permute(vel_node_right, [2 1]), [], 1); % 交错重塑
            vel_node_left_list = reshape(permute(vel_node_left, [2 1]), [], 1); % 交错重塑
            % vel_node_list = reshape(permute(vel_node, [2 1]), [], 1); % 交错重塑
            % vel_node_list=vel_node(:);
        else
            % vel_node_list = [];
            vel_node_right_list = [];
            vel_node_left_list = [];
        end

        % if ~isempty(loc_vel)
        %     vel_node = FindNodewithLocation(loc_vel, node_loc, nodeondeck);
        %     vel_node_list = reshape(permute(vel_node, [2 1]), [], 1);
        % else
        %     vel_node_list = [];
        % end

        if ~isempty(loc_dis)
            dis_node = FindNodewithLocation(loc_dis, node_loc, nodeondeck);
            dis_node_right = dis_node.right;
            dis_node_left = dis_node.left;
            dis_node_weights = dis_node.weights;
            dis_node_right_list = reshape(permute(dis_node_right, [2 1]), [], 1); % 交错重塑
            dis_node_left_list = reshape(permute(dis_node_left, [2 1]), [], 1); % 交错重塑
            % dis_node_list = reshape(permute(dis_node, [2 1]), [], 1); % 交错重塑
            % dis_node_list=dis_node(:);
        else
            % dis_node_list = [];
            dis_node_right_list = [];
            dis_node_left_list = [];
        end

        % if ~isempty(loc_dis)
        %     dis_node = FindNodewithLocation(loc_dis, node_loc, nodeondeck);
        %     dis_node_list = reshape(permute(dis_node, [2 1]), [], 1);
        % else
        %     dis_node_list = [];
        % end




        acc_matrix_left_seq = node2matrixseq(acc_node_left_list, Mapping_data);
        vel_matrix_left_seq = node2matrixseq(vel_node_left_list, Mapping_data);
        dis_matrix_left_seq = node2matrixseq(dis_node_left_list, Mapping_data);
        acc_matrix_right_seq = node2matrixseq(acc_node_right_list, Mapping_data);
        vel_matrix_right_seq = node2matrixseq(vel_node_right_list, Mapping_data);
        dis_matrix_right_seq = node2matrixseq(dis_node_right_list, Mapping_data);

        n_acc = length(acc_node_right_list);
        n_vel = length(vel_node_right_list);
        n_dis = length(dis_node_right_list);
        n_sensors = n_acc + n_vel + n_dis;

        S_a = zeros(n_sensors, size(phi, 1));
        S_v = zeros(n_sensors, size(phi, 1));
        S_d = zeros(n_sensors, size(phi, 1));

        n = n_sensors/(length(loc_dis)+length(loc_vel)+length(loc_acc));

        for k1 = 1:n_sensors
            
            if k1 <= n_acc
                % S_a(k1, acc_matrix_seq(k1)) = 1;
                k2 = ceil(k1/2);
                S_a(k1, acc_matrix_left_seq(k1)) = acc_node_weights(1,k2);
                S_a(k1, acc_matrix_right_seq(k1)) = acc_node_weights(2,k2);
            elseif k1 <= n_acc + n_vel
                k2 = ceil((k1 - n_acc)/2);
                S_v(k1, vel_matrix_left_seq(k1 - n_acc)) = vel_node_weights(1,k2);
                S_v(k1, vel_matrix_right_seq(k1 - n_acc)) = vel_node_weights(2,k2);
            else
                k2 = ceil((k1 - n_acc - n_vel)/2);
                S_d(k1, dis_matrix_left_seq(k1 - n_acc - n_vel)) = dis_node_weights(1,k2);
                S_d(k1, dis_matrix_right_seq(k1 - n_acc - n_vel)) = dis_node_weights(2,k2);
            end

        end

    end
