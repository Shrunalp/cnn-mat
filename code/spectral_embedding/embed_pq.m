function [data] = embed_pq(file_name, n_dim, DREAM3D_order)

% embed_pq will convert Phase Quaternion data from DREAM3D and uses
% get_spectral to find the first n_dim coordinates of the quaternion 
% in the spectral embedding.
% 
% Inputs:
%   file_name     - four components of the normalized quaternion q in the form of a
%                   row vector
%
%   n_dim         - number of coordinates of the quaternion q in the spectral
%                   embedding.
%
%   DREAM3D_order - boolean that indicates the quaternion uses DREAM3D
%                   vector-scalar ordering.
% 
% Outputs:
%   data - first n_dim coordinates of the quaternion q in the spectral
%            embedding in the form of a row vector


    clear all;
    
    PQ = readmatrix(file_name);
    data = zeros(size(PQ,1),n_dim);
    symm = collect_symm(n_dim);
    
    for i = 1:size(PQ,1)
        n = get_spectral(PQ(i,:), n_dim, symm, DREAM3D_order);
        data(i,:) = n;
    end
end



