// Non-manifold discrete surface arrangement and harmonic overlap removal
// using built-in GEO geometry + Plugin(ArrangeSurfaces)

Geometry.OldNewReg = 0;

Macro MakeBox
  p1 = newp; Point(p1) = {xmin, ymin, zmin, lc};
  p2 = newp; Point(p2) = {xmax, ymin, zmin, lc};
  p3 = newp; Point(p3) = {xmax, ymax, zmin, lc};
  p4 = newp; Point(p4) = {xmin, ymax, zmin, lc};
  p5 = newp; Point(p5) = {xmin, ymin, zmax, lc};
  p6 = newp; Point(p6) = {xmax, ymin, zmax, lc};
  p7 = newp; Point(p7) = {xmax, ymax, zmax, lc};
  p8 = newp; Point(p8) = {xmin, ymax, zmax, lc};

  l1 = newc; Line(l1) = {p1, p2};
  l2 = newc; Line(l2) = {p2, p3};
  l3 = newc; Line(l3) = {p3, p4};
  l4 = newc; Line(l4) = {p4, p1};
  l5 = newc; Line(l5) = {p5, p6};
  l6 = newc; Line(l6) = {p6, p7};
  l7 = newc; Line(l7) = {p7, p8};
  l8 = newc; Line(l8) = {p8, p5};
  l9 = newc; Line(l9) = {p1, p5};
  l10 = newc; Line(l10) = {p2, p6};
  l11 = newc; Line(l11) = {p3, p7};
  l12 = newc; Line(l12) = {p4, p8};

  cl1 = newcl; Curve Loop(cl1) = {-l1, -l4, -l3, -l2};
  s1 = news; Plane Surface(s1) = {cl1};
  cl2 = newcl; Curve Loop(cl2) = {l5, l6, l7, l8};
  s2 = news; Plane Surface(s2) = {cl2};
  cl3 = newcl; Curve Loop(cl3) = {l1, l10, -l5, -l9};
  s3 = news; Plane Surface(s3) = {cl3};
  cl4 = newcl; Curve Loop(cl4) = {l2, l11, -l6, -l10};
  s4 = news; Plane Surface(s4) = {cl4};
  cl5 = newcl; Curve Loop(cl5) = {l3, l12, -l7, -l11};
  s5 = news; Plane Surface(s5) = {cl5};
  cl6 = newcl; Curve Loop(cl6) = {l4, l9, -l8, -l12};
  s6 = news; Plane Surface(s6) = {cl6};

  sl = newsl; Surface Loop(sl) = {s1, s2, s3, s4, s5, s6};
  Volume(vTag) = {sl};
Return

// Matrix (Volume 1)
xmin = -1.2; ymin = -1.2; zmin = -1.2;
xmax =  1.2; ymax =  1.2; zmax =  1.2; lc = 0.45; vTag = 1;
Call MakeBox;

// Overlapping equal-priority inclusions (Volumes 2 and 3)
xmin = -0.55; ymin = -0.35; zmin = -0.35;
xmax =  0.20; ymax =  0.35; zmax =  0.35; lc = 0.30; vTag = 2;
Call MakeBox;

xmin = -0.20; ymin = -0.35; zmin = -0.35;
xmax =  0.55; ymax =  0.35; zmax =  0.35; lc = 0.30; vTag = 3;
Call MakeBox;

// RVE clipping window (Volume 4)
xmin = -1.0; ymin = -1.0; zmin = -1.0;
xmax =  1.0; ymax =  1.0; zmax =  1.0; lc = 0.40; vTag = 4;
Call MakeBox;

Physical Volume("Matrix", 1) = {1};
Physical Volume("InclusionA", 2) = {2};
Physical Volume("InclusionB", 3) = {3};

// 1. Mesh closed 2D input shells
Mesh 2;

// 2. Intersect, clip to RVE, partition overlap harmonically, propagate Physicals
Plugin(ArrangeSurfaces).VolumeTags = "1, 2, 3";
Plugin(ArrangeSurfaces).Priorities = "0, 10, 10";
Plugin(ArrangeSurfaces).ClipOutsideVolumeTags = "4";
Plugin(ArrangeSurfaces).ClipPriorities = "20";
Plugin(ArrangeSurfaces).RemoveOverlaps = 1;
Plugin(ArrangeSurfaces).OptimizeSurfaces = 1;
Plugin(ArrangeSurfaces).PropagatePhysicals = 1;
Plugin(ArrangeSurfaces).Run;

// 3. Conforming 3D tetrahedralization of arranged multi-material assembly
Mesh 3;
OptimizeMesh "Gmsh";
