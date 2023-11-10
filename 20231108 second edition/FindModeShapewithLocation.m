%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-11-01 20:54:50
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-11-01 21:01:49
%FilePath: \Exercises-for-Techniques-for-estimation-in-dynamics-systemsf:\git\xihoumen_inverse_force_estimation\20231005 first version\FindModeShapewithLocation.m
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
    mode_deck_loc_two = [Nodeondeck_Info_left.mode_deck,Nodeondeck_Info_right.mode_deck];
    mode_deck_loc = diag(mode_deck_loc_two*weight);
    Modeshape = mode_deck_loc;
end