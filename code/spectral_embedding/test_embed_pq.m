%#ok<*CLALL>
%#ok<*SAGROW>

clear all;

pq_data = "RI_quat_data1024.txt"; %Change file name accordingly
n_dim = 4;

% Iterate through the file to get SEQ for each row vector
for i = 1:length(pq_data)
    currentString = pq_data{i};
    PQ = readmatrix(currentString);

    data = zeros(size(PQ,1),n_dim);
    symm = collect_symm(n_dim);

    for j = 1:size(PQ,1)
        n = get_spectral(PQ(j,2:5), n_dim, symm, true);
        data(j,:) = n;
    end
    file_name = strcat("SE_",currentString):
    writematrix(real(data), file_name);
end
