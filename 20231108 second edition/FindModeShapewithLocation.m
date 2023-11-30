%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-11-01 20:54:50
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-11-30 12:43:14
%FilePath: \Exercises-for-Techniques-for-estimation-in-dynamics-systemsf:\git\xihoumen_inverse_force_estimation\20231108 second edition\FindModeShapewithLocation.m
%Description: 
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Modeshape = FindModeShapewithLocation(loc,node_loc,nodeondeck,KMmapping,nodegap,mode_vec)
    node = FindNodewithLocation(loc,node_loc,nodeondeck);
    node_left = node.left;
    node_right = node.right;
    Nodeondeck_Info_left = Cal_Nodeondeck_Info(node_left,KMmapping,nodegap,mode_vec,'showtext',true);
    Nodeondeck_Info_right = Cal_Nodeondeck_Info(node_right,KMmapping,nodegap,mode_vec,'showtext',true);
    weight = node.weights;

    Nodeondeck_Info_left_mode = Nodeondeck_Info_left.mode_deck;
    Nodeondeck_Info_right_mode = Nodeondeck_Info_right.mode_deck;

    weight_left = weight(1,:);
    weight_right = weight(2,:);

    Nodeondeck_Info_left_mode_addweight = Nodeondeck_Info_left_mode.*weight_left(:);
    Nodeondeck_Info_right_mode_addweight = Nodeondeck_Info_right_mode.*weight_right(:);

    mode_deck_loc = Nodeondeck_Info_left_mode_addweight + Nodeondeck_Info_right_mode_addweight;

    % mode_deck_loc_two = [Nodeondeck_Info_left.mode_deck(:),Nodeondeck_Info_right.mode_deck(:)];
    % mode_deck_loc = mode_deck_loc_two*weight;
    Modeshape = mode_deck_loc;
end