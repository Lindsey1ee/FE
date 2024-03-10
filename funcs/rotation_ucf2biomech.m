function xyz1 = rotation_ucf2biomech(ucf_xyz)

% UCF BRaIN Lab axes
% +x, right; -x, left
% +y, up; -y, down
% +z, posterior; -z, anterior
%
% more intitutive biomech
% +x, right; -x, left
% +y, anterior; -y, posterior
% +z, up; -z, down
%
% rotate 90 degrees about the x-axis

theta = 90;

% 3d rotation about x-axis
rot_x = [1 0 0; ...
    0 cosd(theta) -sind(theta); ...
    0 sind(theta) cosd(theta)];

xyz1 = rot_x*ucf_xyz';
xyz1 = xyz1';
