%#ok<*CLALL>
%#ok<*SAGROW>

clear all;

n_pts = 1000;
n_dim = 12;
n_neigh = 1000;

generators = {angles_to_CK(pi / 2., 0., 0.), ...
              angles_to_CK(pi / 2., pi / 2., 0.)};

TOL = sqrt(eps);

% Generate points on the 3-sphere using the Fibonacci sphere algorithm
n_pts1 = n_pts + 1;
ph = zeros(n_pts1, 1);
th = zeros(n_pts1, 1);
w = zeros(n_pts1, 1);

c1 = 2. * pi;
c2 = pi / n_pts;
for a = 0:n_pts
    ph(a + 1) = c1 * mod(a * sqrt(37.), 1.);
    th(a + 1) = acos(1. - 2. * mod(a * sqrt(421.), 1.));
    w(a + 1) = newton(c2 * a);
end

% Construct coefficients
symm = collect_symm(n_dim);

% Generate symmetries
len0 = 0;
while true
    len1 = length(generators);
    for a = 1:len1
        for b = max(a, len0 + 1):len1
            CK = generators{a} * generators{b};
            include = true;
            for c = 1:length(generators)
                if norm(CK - generators{c}) < TOL || ...
                   norm(CK + generators{c}) < TOL
                    include = false;
                    break;
                end
            end
            if include
                generators{end + 1} = CK;
            end
        end
    end
    if length(generators) == len1
        break;
    else
        len0 = len1;
    end
end
n_gen = length(generators);

% Generate equivalent quaternion sets and spectral embedding
Q = zeros(n_pts1, 4, 2. * n_gen);
S = zeros(n_pts1, n_dim);

for a = 1:n_pts1
    CK = angles_to_CK(w(a), th(a), ph(a));
    for b = 1:n_gen
        [w_s, th_s, ph_s] = CK_to_angles(CK * generators{b});
        q = [cos(w_s / 2.),...
             sin(w_s / 2.) * sin(th_s) * cos(ph_s),...
             sin(w_s / 2.) * sin(th_s) * sin(ph_s),...
             sin(w_s / 2.) * cos(th_s)];
        Q(a, :, b) = q;
        Q(a, :, b + n_gen) = -q;
    end
    S(a, :) = get_spectral(Q(a, :, 1), n_dim, symm, false);
end

% Find spectral embedding distance for points nearby as quaternions
dist_q = zeros(n_neigh * n_pts1, 1);
dist_s = zeros(n_neigh * n_pts1, n_dim);

d = zeros(n_pts1, 1);
for a = 0:n_pts
    rows = (n_neigh * a + 1):(n_neigh * (a + 1));

    Qa = repmat(Q(a + 1, :, 1), 1, 1, 2. * n_gen);
    for b = 0:n_pts
        sft = Q(b + 1, :, :) - Qa;
        d(b + 1) = min(squeeze(sum(sft.^2, 2)));
    end
    [~, I] = mink(d, n_neigh + 1, 'comparisonmethod', 'real');
    dist_q(rows) = d(I(2:(n_neigh + 1)));

    tmp = S(I(2:(n_neigh + 1)), :) - repmat(S(a + 1, :), n_neigh, 1);
    for b = 1:n_dim
        dist_s(rows, b) = sum(tmp(:, 1:b).^2, 2);
    end
end

% Plot the distances as truncated spectral embeddings vs distances as
% quaternions
for a = 1:n_dim
    figure(a), clf;
    plot(dist_q, dist_s(:, a), '.', 'markersize', 10);
    xlabel('Quaternion distance');
    ylabel('Spectral distance');
    title([num2str(a), ' dimensions']);
end

% % Find quaternion distance for points nearby in the spectral embedding
% dist_q = zeros(n_neigh * n_pts1, n_dim);
% dist_s = zeros(n_neigh * n_pts1, n_dim);
% 
% d = zeros(n_pts1, 1);
% for a = 0:n_pts
%     rows = (n_neigh * a + 1):(n_neigh * (a + 1));
% 
%     Qa = repmat(Q(a + 1, :, 1), 1, 1, 2. * n_gen);
%     for b = 1:n_dim
%         tmp = S(:, 1:b) - repmat(S(a + 1, 1:b), n_pts1, 1);
%         d = sum(tmp.^2, 2);
%         [~, I] = mink(d, n_neigh + 1, 'comparisonmethod', 'real');
%         dist_s(rows, b) = d(I(2:(n_neigh + 1)));
% 
%         for c = 2:(n_neigh + 1)
%             tmp = Q(I(c), :, :) - Qa;
%             d(c) = min(squeeze(sum(tmp.^2, 2)));
%         end
%         dist_q(rows, b) = d(2:(n_neigh + 1));
%     end
% end
% 
% % Plot the distances as truncated spectral embeddings vs distances as
% % quaternions
% for a = 1:n_dim
%     figure(a), clf;
%     plot(dist_s(:, a), dist_q(:, a), '.', 'markersize', 10);
%     xlabel('Spectral distance');
%     ylabel('Quaternion distance');
%     title([num2str(a), ' dimensions']);
% end

save('cubic_test.mat', 'w', 'th', 'ph', 'generators', 'dist_q', 'dist_s');
