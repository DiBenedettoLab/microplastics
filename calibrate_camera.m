function rectify_quad = calibrate_camera(I,W)
% undistortion and rectification of points in an image (based on
% Fujita et al 1998, but with a quadratic transformation rather than
% linear)
%
% image2world_coords: function handle to map image coordinates to world coordinates 
% I: set of known calibration points in image coordinates (n x 2 vector) [px]
% W: set of known calibration points in world coordinates (n x 2 vector) [m]

A = [I.^2, I, ones(size(I,1),1), -I(:,1).^2.*W(:,1), -I(:,2).^2.*W(:,1), -I(:,1).*W(:,1), -I(:,2).*W(:,1), zeros(size(I,1),5);
    zeros(size(I,1),5), -I(:,1).^2.*W(:,2), -I(:,2).^2.*W(:,2), -I(:,1).*W(:,2), -I(:,2).*W(:,2), I.^2, I, ones(size(I,1),1)];
Z = [W(:,1); W(:,2)];
B = (A'*A)^-1*A'*Z;

% function to map image coords to world coords
rectify_quad = @(I) [(B(1)*I(:,1).^2 + B(2)*I(:,2).^2 + B(3)*I(:,1) + B(4)*I(:,2) + B(5))./ ...
    (B(6)*I(:,1).^2 + B(7)*I(:,2).^2 + B(8)*I(:,1) + B(9)*I(:,2) + 1), ...
    (B(10)*I(:,1).^2 + B(11)*I(:,2).^2 + B(12)*I(:,1) + B(13)*I(:,2) + B(14))./ ...
    (B(6)*I(:,1).^2 + B(7)*I(:,2).^2 + B(8)*I(:,1) + B(9)*I(:,2) + 1)];

end