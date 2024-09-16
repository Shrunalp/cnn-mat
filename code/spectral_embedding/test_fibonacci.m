% This is a test of the Fibonacci sphere algorithm for distributing points on
% the 2-sphere and 3-sphere. Shows both the distribution of points and the 
% histogram of the Euclidean distances in the 10-nearest neighbor graphs.

N = 10000;
N1 = N + 1;
p = zeros(N1, 1);
t = zeros(N1, 1);
o = zeros(N1, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize lattice on 2-sphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gold = (1. + sqrt(5.)) / 2.;
% 
% c1 = 2. * pi;
% c2 = 2. / N;
% for a = 0:N
%     p(a + 1) = c1 * mod(a * gold, 1.);
%     t(a + 1) = acos(1. - c2 * a);
% end
% 
% pts = [sin(t) .* cos(p), ...
%        sin(t) .* sin(p), ...
%        cos(t)];
% 
% dist = zeros(10 * N1, 1);
% for a = 0:N
%     sft = pts - repmat(pts(a + 1, :), N1, 1);
%     d = sort(sum(sft.^2, 2));
%     dist((10 * a + 1):(10 * (a + 1))) = d(2:11);
% end
% 
% figure(1), clf;
% [y, x] = ecdf(dist);
% plot(x, y);
% 
% [X, Y, Z] = sphere(100);
% figure(2), clf;
% hold on, box on, grid on;
% surf(X, Y, Z, 'edgealpha', 0, 'facecolor', [0., 0., 1.], 'facealpha', 0.1);
% plot3(pts(:, 1), pts(:, 2), pts(:, 3), '.k', 'markersize', 10.);
% axis equal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Visualize lattice on 3-sphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c1 = 2. * pi;
c2 = pi / N;
for a = 0:N
    p(a + 1) = c1 * mod(a * sqrt(37.), 1.);
    t(a + 1) = acos(1. - 2. * mod(a * sqrt(421.), 1.));
    o(a + 1) = newton(c2 * a);
end

pts = [cos(o), ...
       sin(o) .* sin(t) .* cos(p), ...
       sin(o) .* sin(t) .* sin(p), ...
       sin(o) .* cos(t)];

dist = zeros(10 * N1, 1);
for a = 0:N
    sft = pts - repmat(pts(a + 1, :), N1, 1);
    sft = sft(all(abs(sft) < 0.2, 2), :);
    d = sort(sum(sft.^2, 2));
    dist((10 * a + 1):(10 * (a + 1))) = d(2:11);
end

figure(1), clf;
[y, x] = ecdf(dist);
plot(x, y);

[X, Y, Z] = sphere(100);
figure(2), clf;
hold on, box on, grid on;
surf(X, Y, Z, 'edgealpha', 0, 'facecolor', [0., 0., 1.], 'facealpha', 0.1);
f = (pts(:, 1) < 0.02) & (pts(:, 1) > 0.0);
plot3(pts(f, 2), pts(f, 3), pts(f, 4), '.k', 'markersize', 10.);
axis equal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate different irrational numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,...
%           59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, ...
%           127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, ...
%           191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, ...
%           257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, ...
%           331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, ...
%           401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, ...
%           467, 479, 487, 491, 499, 503, 509, 521, 523, 541];
% primes = sqrt(primes);
% np = numel(primes);
% 
% obj1 = zeros(np, np);
% obj2 = zeros(np, np);
% obj3 = zeros(np, np);
% c2 = pi / N;
% for a = 0:N
%     o(a + 1) = newton(c2 * a);
% end
%         
% for z1 = 1:np
%     for z2 = [1:(z1 - 1), (z1 + 1):np]
%         c1 = 2. * pi;
%         for a = 0:N
%             p(a + 1) = c1 * mod(a * primes(z1), 1.);
%             t(a + 1) = acos(1. - 2. * mod(a * primes(z2), 1.));
%         end
% 
%         pts = [cos(o), ...
%                sin(o) .* sin(t) .* cos(p), ...
%                sin(o) .* sin(t) .* sin(p), ...
%                sin(o) .* cos(t)];
% 
%         dist = zeros(10 * N1, 1);
%         for a = 0:N
%             sft = pts - repmat(pts(a + 1, :), N1, 1);
%             sft = sft(all(abs(sft) < 0.2, 2), :);
%             d = sort(sum(sft.^2, 2));
%             dist((10 * a + 1):(10 * (a + 1))) = d(2:11);
%         end
% 
%         dist = sort(dist);
%         obj1(z1, z2) = dist(100);
%         obj2(z1, z2) = dist(1000);
%         obj3(z1, z2) = dist(10000);
%     end
% end
