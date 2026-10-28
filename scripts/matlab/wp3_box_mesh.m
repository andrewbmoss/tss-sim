function TR = wp3_box_mesh(bx, h)
%WP3_BOX_MESH  Closed triangulated box [-bx1/2,bx1/2]x[-bx2/2,bx2/2]x[0,bx3], edge about h.
n = max(2, ceil(bx ./ h)); P = []; T = [];
faces = {[2 3 1 0], [2 3 1 1], [1 3 2 0], [1 3 2 1], [1 2 3 0], [1 2 3 1]};
lo = [-bx(1)/2, -bx(2)/2, 0]; hi = [bx(1)/2, bx(2)/2, bx(3)];
for k = 1:6
  a = faces{k}(1); b = faces{k}(2); c = faces{k}(3); s = faces{k}(4);
  [u, v] = ndgrid(linspace(lo(a), hi(a), n(a)+1), linspace(lo(b), hi(b), n(b)+1));
  q = zeros(numel(u), 3); q(:, a) = u(:); q(:, b) = v(:); q(:, c) = lo(c) + s*(hi(c)-lo(c));
  na = n(a)+1; nb = n(b)+1; base = size(P, 1); [i, j] = ndgrid(1:n(a), 1:n(b));
  i1 = sub2ind([na nb], i(:), j(:)); i2 = sub2ind([na nb], i(:)+1, j(:));
  i3 = sub2ind([na nb], i(:)+1, j(:)+1); i4 = sub2ind([na nb], i(:), j(:)+1);
  T = [T; base + [i1 i2 i3]; base + [i1 i3 i4]]; P = [P; q];
end
[P, ~, ic] = uniquetol(P, 1e-9, 'ByRows', true); T = ic(T);
c0 = (lo + hi) / 2; nrm = cross(P(T(:,2),:) - P(T(:,1),:), P(T(:,3),:) - P(T(:,1),:), 2);
cen = (P(T(:,1),:) + P(T(:,2),:) + P(T(:,3),:)) / 3; flip = sum(nrm .* (cen - c0), 2) < 0;
T(flip, :) = T(flip, [1 3 2]);  % outward normals
TR = triangulation(T, P);
end
