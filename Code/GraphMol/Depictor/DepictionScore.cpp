//
//  Experimental high-quality 2D depiction support.
//
#include "DepictionScore.h"

#include <GraphMol/ROMol.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Conformer.h>
#include <Geometry/point.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace RDDepict {
namespace {

constexpr double kEps = 1.0e-8;
constexpr double kPi = 3.14159265358979323846;
constexpr double kAtomClashDistance = 0.55;
constexpr double kAtomBondClashDistance = 0.28;
constexpr double kCongestionDistance = 1.05;
constexpr double kMinPreferredAngleDeg = 28.0;

RDGeom::Point2D point2D(const RDKit::Conformer &conf, unsigned int idx) {
  const auto &p = conf.getAtomPos(idx);
  return RDGeom::Point2D(p.x, p.y);
}

double cross(const RDGeom::Point2D &a, const RDGeom::Point2D &b,
             const RDGeom::Point2D &c) {
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

bool strictSegmentIntersection(const RDGeom::Point2D &a,
                               const RDGeom::Point2D &b,
                               const RDGeom::Point2D &c,
                               const RDGeom::Point2D &d) {
  const auto c1 = cross(a, b, c);
  const auto c2 = cross(a, b, d);
  const auto c3 = cross(c, d, a);
  const auto c4 = cross(c, d, b);
  return ((c1 > kEps && c2 < -kEps) || (c1 < -kEps && c2 > kEps)) &&
         ((c3 > kEps && c4 < -kEps) || (c3 < -kEps && c4 > kEps));
}

double pointSegmentDistance(const RDGeom::Point2D &p,
                            const RDGeom::Point2D &a,
                            const RDGeom::Point2D &b) {
  const double vx = b.x - a.x;
  const double vy = b.y - a.y;
  const double wx = p.x - a.x;
  const double wy = p.y - a.y;
  const double vv = vx * vx + vy * vy;
  if (vv < kEps) {
    return std::hypot(wx, wy);
  }
  double t = (wx * vx + wy * vy) / vv;
  t = std::clamp(t, 0.0, 1.0);
  const double dx = p.x - (a.x + t * vx);
  const double dy = p.y - (a.y + t * vy);
  return std::hypot(dx, dy);
}

double angleDeg(const RDGeom::Point2D &a, const RDGeom::Point2D &center,
                const RDGeom::Point2D &b) {
  const double ax = a.x - center.x;
  const double ay = a.y - center.y;
  const double bx = b.x - center.x;
  const double by = b.y - center.y;
  const double al = std::hypot(ax, ay);
  const double bl = std::hypot(bx, by);
  if (al < kEps || bl < kEps) {
    return 0.0;
  }
  double cosine = (ax * bx + ay * by) / (al * bl);
  cosine = std::clamp(cosine, -1.0, 1.0);
  return std::acos(cosine) * 180.0 / kPi;
}

}  // namespace

DepictionScore scoreDepiction(const RDKit::ROMol &mol, int confId,
                              const DepictionScoreWeights &weights) {
  DepictionScore result;
  if (!mol.getNumAtoms() || !mol.getNumConformers()) {
    result.total = std::numeric_limits<double>::infinity();
    return result;
  }

  const auto &conf = mol.getConformer(confId);

  // Materialize the C++20 bond range once. This avoids the legacy
  // ConstBondIterator API, whose definition is intentionally not pulled in by
  // ROMol.h, and also makes the O(E^2) crossing scan straightforward.
  std::vector<const RDKit::Bond *> bonds;
  bonds.reserve(mol.getNumBonds(false));
  for (const auto *bond : mol.bonds()) {
    bonds.push_back(bond);
  }

  for (size_t i = 0; i < bonds.size(); ++i) {
    const auto *b1 = bonds[i];
    for (size_t j = i + 1; j < bonds.size(); ++j) {
      const auto *b2 = bonds[j];
      const auto a1 = b1->getBeginAtomIdx();
      const auto a2 = b1->getEndAtomIdx();
      const auto c1 = b2->getBeginAtomIdx();
      const auto c2 = b2->getEndAtomIdx();
      if (a1 == c1 || a1 == c2 || a2 == c1 || a2 == c2) {
        continue;
      }
      if (strictSegmentIntersection(point2D(conf, a1), point2D(conf, a2),
                                    point2D(conf, c1), point2D(conf, c2))) {
        ++result.bondCrossings;
      }
    }
  }

  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    const auto pi = point2D(conf, i);
    for (unsigned int j = i + 1; j < mol.getNumAtoms(); ++j) {
      if (mol.getBondBetweenAtoms(i, j)) {
        continue;
      }
      const auto pj = point2D(conf, j);
      const double d = std::hypot(pi.x - pj.x, pi.y - pj.y);
      if (d < kAtomClashDistance) {
        ++result.atomAtomClashes;
      }
      if (d < kCongestionDistance && d > kEps) {
        const double x = (kCongestionDistance - d) / kCongestionDistance;
        result.congestionPenalty += x * x;
      }
    }
  }

  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    const auto p = point2D(conf, i);
    for (const auto *bond : bonds) {
      const auto a = bond->getBeginAtomIdx();
      const auto b = bond->getEndAtomIdx();
      if (i == a || i == b || mol.getBondBetweenAtoms(i, a) ||
          mol.getBondBetweenAtoms(i, b)) {
        continue;
      }
      if (pointSegmentDistance(p, point2D(conf, a), point2D(conf, b)) <
          kAtomBondClashDistance) {
        ++result.atomBondClashes;
      }
    }
  }

  for (const auto *atom : mol.atoms()) {
    std::vector<unsigned int> nbrs;
    for (const auto *nbr : mol.atomNeighbors(atom)) {
      nbrs.push_back(nbr->getIdx());
    }
    if (nbrs.size() < 2) {
      continue;
    }
    const auto center = point2D(conf, atom->getIdx());
    for (size_t i = 0; i < nbrs.size(); ++i) {
      for (size_t j = i + 1; j < nbrs.size(); ++j) {
        const double angle = angleDeg(point2D(conf, nbrs[i]), center,
                                      point2D(conf, nbrs[j]));
        if (angle < kMinPreferredAngleDeg) {
          const double x = (kMinPreferredAngleDeg - angle) /
                           kMinPreferredAngleDeg;
          result.anglePenalty += x * x;
        }
      }
    }
  }

  std::vector<double> bondLengths;
  bondLengths.reserve(bonds.size());
  for (const auto *bond : bonds) {
    const auto p1 = point2D(conf, bond->getBeginAtomIdx());
    const auto p2 = point2D(conf, bond->getEndAtomIdx());
    bondLengths.push_back(std::hypot(p1.x - p2.x, p1.y - p2.y));
  }
  if (!bondLengths.empty()) {
    const double mean = std::accumulate(bondLengths.begin(), bondLengths.end(),
                                        0.0) /
                        bondLengths.size();
    if (mean > kEps) {
      for (double length : bondLengths) {
        const double rel = (length - mean) / mean;
        result.bondLengthPenalty += rel * rel;
      }
      result.bondLengthPenalty /= bondLengths.size();
    }
  }

  double minX = std::numeric_limits<double>::max();
  double maxX = std::numeric_limits<double>::lowest();
  double minY = std::numeric_limits<double>::max();
  double maxY = std::numeric_limits<double>::lowest();
  for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    const auto p = point2D(conf, i);
    minX = std::min(minX, p.x);
    maxX = std::max(maxX, p.x);
    minY = std::min(minY, p.y);
    maxY = std::max(maxY, p.y);
  }
  const double width = std::max(maxX - minX, kEps);
  const double height = std::max(maxY - minY, kEps);
  const double ratio = std::max(width / height, height / width);
  if (ratio > 3.0) {
    result.aspectRatioPenalty = std::log1p(ratio - 3.0);
  }

  result.total = weights.bondCrossing * result.bondCrossings +
                 weights.atomAtomClash * result.atomAtomClashes +
                 weights.atomBondClash * result.atomBondClashes +
                 weights.congestion * result.congestionPenalty +
                 weights.angle * result.anglePenalty +
                 weights.bondLength * result.bondLengthPenalty +
                 weights.aspectRatio * result.aspectRatioPenalty;
  return result;
}

}  // namespace RDDepict
