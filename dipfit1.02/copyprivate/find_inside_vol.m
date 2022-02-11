function [inside, outside] = find_inside_vol(pos, vol);

% FIND_INSIDE_VOL locates dipole locations inside/outside the brain
% compartment of a volume conductor model.
% 
% [inside, outside] = find_inside_vol(pos, vol)
%
% where
%   pos      Nx3 matrix with dipole positions
%   vol      structure with volume conductor model
%   inside   index of all dipoles inside the brain compartment
%   outside  index of all dipoles outside the brain compartment

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: find_inside_vol.m,v $
% Revision 1.3  2003/08/04 09:19:30  roberto
% fixed bug for dipole on BEM volume boundary
%
% Revision 1.2  2003/03/24 12:30:06  roberto
% added support for multi-sphere volume conductor model
%

% single-sphere volume conductor model
if isfield(vol, 'r') & all(size(vol.r)==[1 1])
  if isfield(vol, 'o')
    % temporary shift dipole positions toward origin
    tmp = pos - repmat(vol.o, size(pos,1), 1);
  else
    tmp = pos;
  end
  outside = find(sqrt(sum(tmp.^2, 2))>min(vol.r));
  inside  = setdiff(1:size(pos,1),outside)';

% multi-sphere volume conductor model
elseif isfield(vol, 'r') & ~all(size(vol.r)==[1 1])

%  nspheres = size(vol.r,1);
%  ndipoles = size(pos,1);
%  inside = zeros(ndipoles,1);
%  for sph=1:nspheres
%  for dip=1:ndipoles
%    if inside(dip)
%      % the dipole has already been detected in one of the other spheres
%      continue
%    end
%    inside(dip) = (norm(pos(dip,:) - vol.o(sph,:)) <= vol.r(sph));
%  end
%  end
%  outside = find(inside==0);
%  inside  = find(inside==1);

  % this is a much faster implementation
  nspheres = size(vol.r,1);
  ndipoles = size(pos,1);
  inside = zeros(ndipoles,1);
  for sph=1:nspheres
    % temporary shift dipole positions toward origin
    if isfield(vol, 'o')
      tmp = pos - repmat(vol.o(sph,:), [ndipoles 1]);
    else
      tmp = pos;
    end
    flag = (sqrt(sum(tmp.^2,2)) <= vol.r(sph));
    inside = inside + flag;
  end
  outside = find(inside==0);
  inside  = find(inside>0);

% realistic BEM volume conductor model
elseif isfield(vol, 'bnd')
  if isfield(vol, 'source')
    % use the specified source compartment
    pnt = vol.bnd(vol.source).pnt;
    tri = vol.bnd(vol.source).tri;
  else
    % locate the innermost compartment and remember it
    vol.source = find_innermost_boundary(vol.bnd)
    pnt = vol.bnd(vol.source).pnt;
    tri = vol.bnd(vol.source).tri;
  end
  % determine the dipole positions that are inside the brain compartment
  inside  = find(bounding_mesh(pos, pnt, tri)==1)';
  outside = setdiff(1:size(pos,1),inside)';

% unrecognized volume conductor model
else
  error('unrecognized volume conductor model');
end

% ensure column vectors
inside = inside(:);
outside = outside(:);

