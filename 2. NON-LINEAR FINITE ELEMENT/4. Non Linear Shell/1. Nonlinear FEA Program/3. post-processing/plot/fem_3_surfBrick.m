function [qu] = fem_3_surfBrick(p,t)

% -------------------------------------------------------------------------

FacesA = [t(:,[1,2,6,5]);
          t(:,[1,5,8,4]);
          t(:,[5,6,7,8]);
          t(:,[2,6,7,3]);
          t(:,[1,2,3,4]);
          t(:,[4,3,7,8])];

Faces = sort(FacesA,2);
[foo,ix,jx] = unique(Faces,'rows');
vec   = histc(jx,1:max(jx));
qx    = find(vec==1);
qu    = FacesA(ix(qx),:);
