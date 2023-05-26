%% Funkcja generująca dane dla kinematyki w przód dla pary obrotowej
function [r_0_i, dr_0_i ,ddr_0_i ,R_0_i, w_0_i, dw_0_i, r_com_i, dr_com_i, ddr_com_i, u] = forward_kinematics_obr(r_0_i_1,dr_0_i_1,ddr_0_i_1,R_0_i_1,w_0_i_1,dw_0_i_1,s_com,type,q,dq,ddq,r_const,R_const)

if type == "z"
    w = [0,0,1]';
    [r_0_i, dr_0_i ,ddr_0_i ,R_0_i, w_0_i, dw_0_i, r_com_i, dr_com_i, ddr_com_i, u] = forward_kinematics(r_0_i_1,dr_0_i_1,ddr_0_i_1,R_0_i_1,w_0_i_1,dw_0_i_1,r_const,zeros(3,1),zeros(3,1),R_const*RotZ(q),R_const*w*dq,R_const*w*ddq,s_com, w);
end

if type == "y"
    w = [0,1,0]';
    [r_0_i, dr_0_i ,ddr_0_i ,R_0_i, w_0_i, dw_0_i, r_com_i, dr_com_i, ddr_com_i, u] = forward_kinematics(r_0_i_1,dr_0_i_1,ddr_0_i_1,R_0_i_1,w_0_i_1,dw_0_i_1,r_const,zeros(3,1),zeros(3,1),R_const*RotY(q),R_const*w*dq,R_const*w*ddq,s_com, w);
end

if type == "x"
    w = [1,0,0]';
    [r_0_i, dr_0_i ,ddr_0_i ,R_0_i, w_0_i, dw_0_i, r_com_i, dr_com_i, ddr_com_i, u] = forward_kinematics(r_0_i_1,dr_0_i_1,ddr_0_i_1,R_0_i_1,w_0_i_1,dw_0_i_1,r_const,zeros(3,1),zeros(3,1),R_const*RotX(q),R_const*w*dq,R_const*w*ddq,s_com, w);
end

end

