function [T] = point_registration(P_source,P_end)

% Pa, Pb must be homogeneous coordinates with dimension 3/4 x n, where n is
% the number of points

P_source = P_source(1:3,:); % Remove 1s
% P_source(4,:)=[] this is an alternative way to remove a line
P_end = P_end(1:3,:);

dims = size(P_source);
n_points = dims(2);  % Number of points

B_source = mean(P_source,2);    % Center of mass
B_end = mean(P_end,2);

% Cross correlation function
cross_matrix = zeros(dims(1));
for i = 1 : n_points
    cross_matrix = cross_matrix + ((P_source(:,i) - B_source) * (P_end(:,i) - B_end)');
end
cross_matrix = 1 / n_points * cross_matrix;

% Calculate the parameters for the quaternion extraction
A = cross_matrix - cross_matrix';
D = [A(2,3) A(3,1) A(1,2)]';
Q = [[trace(cross_matrix), D']
    [D, (cross_matrix + cross_matrix' - trace(cross_matrix) * eye(3))]];

% Eigenval/Eigenvect extraction
[E_vec,E_val] = eig(Q);

% Find the maximum eigenval and its position
[max_v max_p] = max(diag(E_val));


% extract the eigenvector corresponding to max_p which will be our
% quaternion. Shape should be [1,4]
quat = E_vec(:,max_p)';

% NOT NECESSARY: convert to quaternion format
%quat = Quaternion(quat);
%quat = quat.double;

% now we need the rotation matrix
% Extract the rotation matrix
T = eye(4);
T(1:3,1:3) = quat2rotm(quat);

% Calculates the translation component according to eq. 26 (check the
% shapes!)
T(1:3,4) = B_end - (T(1:3,1:3)) * B_source;
