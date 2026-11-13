#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>
#include <gmsh.h>

static int addGeoBox(double xmin, double ymin, double zmin, double xmax,
                     double ymax, double zmax, double lc)
{
  int p[8];
  p[0] = gmsh::model::geo::addPoint(xmin, ymin, zmin, lc);
  p[1] = gmsh::model::geo::addPoint(xmax, ymin, zmin, lc);
  p[2] = gmsh::model::geo::addPoint(xmax, ymax, zmin, lc);
  p[3] = gmsh::model::geo::addPoint(xmin, ymax, zmin, lc);
  p[4] = gmsh::model::geo::addPoint(xmin, ymin, zmax, lc);
  p[5] = gmsh::model::geo::addPoint(xmax, ymin, zmax, lc);
  p[6] = gmsh::model::geo::addPoint(xmax, ymax, zmax, lc);
  p[7] = gmsh::model::geo::addPoint(xmin, ymax, zmax, lc);

  int l[12];
  l[0] = gmsh::model::geo::addLine(p[0], p[1]);
  l[1] = gmsh::model::geo::addLine(p[1], p[2]);
  l[2] = gmsh::model::geo::addLine(p[2], p[3]);
  l[3] = gmsh::model::geo::addLine(p[3], p[0]);
  l[4] = gmsh::model::geo::addLine(p[4], p[5]);
  l[5] = gmsh::model::geo::addLine(p[5], p[6]);
  l[6] = gmsh::model::geo::addLine(p[6], p[7]);
  l[7] = gmsh::model::geo::addLine(p[7], p[4]);
  l[8] = gmsh::model::geo::addLine(p[0], p[4]);
  l[9] = gmsh::model::geo::addLine(p[1], p[5]);
  l[10] = gmsh::model::geo::addLine(p[2], p[6]);
  l[11] = gmsh::model::geo::addLine(p[3], p[7]);

  const int loops[6][4] = {
    {-l[0], -l[3], -l[2], -l[1]}, {l[4], l[5], l[6], l[7]},
    {l[0], l[9], -l[4], -l[8]},   {l[1], l[10], -l[5], -l[9]},
    {l[2], l[11], -l[6], -l[10]}, {l[3], l[8], -l[7], -l[11]}};

  std::vector<int> faces;
  for(int i = 0; i < 6; ++i) {
    std::vector<int> lp(loops[i], loops[i] + 4);
    int cl = gmsh::model::geo::addCurveLoop(lp);
    faces.push_back(gmsh::model::geo::addPlaneSurface(std::vector<int>(1, cl)));
  }
  int sl = gmsh::model::geo::addSurfaceLoop(faces);
  return gmsh::model::geo::addVolume(std::vector<int>(1, sl));
}

static void checkVolumeAndQuality(const std::vector<int> &volTags,
                                  double expectedVol)
{
  std::vector<std::size_t> nodeTags;
  std::vector<double> nodeCoords, paramCoords;
  gmsh::model::mesh::getNodes(nodeTags, nodeCoords, paramCoords);
  std::map<std::size_t, std::size_t> tagIdx;
  for(std::size_t i = 0; i < nodeTags.size(); ++i) tagIdx[nodeTags[i]] = i;

  double totalVol = 0.0;
  double minSicn = 1.0;
  for(std::size_t v = 0; v < volTags.size(); ++v) {
    std::vector<int> etypes;
    std::vector<std::vector<std::size_t> > etags, enodes;
    gmsh::model::mesh::getElements(etypes, etags, enodes, 3, volTags[v]);
    for(std::size_t k = 0; k < etypes.size(); ++k) {
      if(etypes[k] != 4) continue;
      std::vector<double> quals;
      gmsh::model::mesh::getElementQualities(etags[k], quals, "minSICN");
      for(std::size_t q = 0; q < quals.size(); ++q)
        minSicn = std::min(minSicn, quals[q]);
      for(std::size_t j = 0; j < enodes[k].size(); j += 4) {
        const double *p0 = &nodeCoords[3 * tagIdx[enodes[k][j]]];
        const double *p1 = &nodeCoords[3 * tagIdx[enodes[k][j + 1]]];
        const double *p2 = &nodeCoords[3 * tagIdx[enodes[k][j + 2]]];
        const double *p3 = &nodeCoords[3 * tagIdx[enodes[k][j + 3]]];
        const double ux = p1[0] - p0[0], uy = p1[1] - p0[1], uz = p1[2] - p0[2];
        const double vx = p2[0] - p0[0], vy = p2[1] - p0[1], vz = p2[2] - p0[2];
        const double wx = p3[0] - p0[0], wy = p3[1] - p0[1], wz = p3[2] - p0[2];
        const double det = ux * (vy * wz - vz * wy) -
                           uy * (vx * wz - vz * wx) +
                           uz * (vx * wy - vy * wx);
        totalVol += std::fabs(det) / 6.0;
      }
    }
  }
  if(std::fabs(totalVol - expectedVol) > 1e-5 * expectedVol || minSicn <= 0.0) {
    std::fprintf(stderr,
                 "Verification failed: vol=%.8f (expected %.8f), minSICN=%g\n",
                 totalVol, expectedVol, minSicn);
    std::exit(1);
  }
}

int main(int argc, char **argv)
{
  gmsh::initialize(argc, argv);
  gmsh::option::setNumber("General.Verbosity", 2);

  gmsh::model::add("arrange_harmonic_cpp");
  int v0 = addGeoBox(-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 0.45);
  int v1 = addGeoBox(-0.55, -0.35, -0.40, 0.35, 0.35, 0.40, 0.32);
  int v2 = addGeoBox(-0.20, -0.55, -0.35, 0.55, 0.25, 0.45, 0.32);
  int v3 = addGeoBox(-0.30, -0.25, -0.55, 0.30, 0.55, 0.55, 0.32);
  gmsh::model::geo::synchronize();

  gmsh::model::mesh::generate(2);

  gmsh::vectorpair objects;
  objects.push_back(std::make_pair(3, v0));
  objects.push_back(std::make_pair(3, v1));
  objects.push_back(std::make_pair(3, v2));
  objects.push_back(std::make_pair(3, v3));

  std::vector<int> prios;
  prios.push_back(0);
  prios.push_back(10);
  prios.push_back(10);
  prios.push_back(10);

  gmsh::vectorpair outVols;
  std::vector<gmsh::vectorpair> outMap;
  gmsh::model::mesh::arrangeSurfaces(objects, gmsh::vectorpair(), outVols,
                                     outMap, prios, std::vector<int>(),
                                     std::vector<int>(), true, true);

  std::vector<int> allChildren;
  const char *names[4] = {"Matrix", "PhaseA", "PhaseB", "PhaseC"};
  for(std::size_t i = 0; i < 4; ++i) {
    std::vector<int> ctags;
    for(std::size_t k = 0; k < outMap[i].size(); ++k) {
      if(outMap[i][k].first == 3) {
        ctags.push_back(outMap[i][k].second);
        allChildren.push_back(outMap[i][k].second);
      }
    }
    if(ctags.empty()) std::exit(1);
    gmsh::model::addPhysicalGroup(3, ctags, static_cast<int>(i + 1), names[i]);
  }

  gmsh::model::mesh::generate(3);
  gmsh::model::mesh::optimize("", true);
  checkVolumeAndQuality(allChildren, 8.0);

  bool popup = true;
  for(int i = 1; i < argc; ++i) {
    if(std::strcmp(argv[i], "-nopopup") == 0) popup = false;
  }
  if(popup) gmsh::fltk::run();

  gmsh::finalize();
  return 0;
}
