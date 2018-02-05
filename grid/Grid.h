#ifndef GRID_H
#define GRID_H

/*
* GridMC - Grid based MC simulation package.
*
* Copyright 2012, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <vector>
#include <utility>
#include "../util/global.h"
#include "../util/IntVector.h"
#include "../util/Vector.h"
#include "PoissonGreenFunc.h"

namespace GridMC
{

   using namespace Util;

   /**
   * An Grid object represents an empty grid.
   *
   * \ingroup Grid_Module
   */
   class Grid 
   {
   
   public:
  
      /**
      * Default constructor.
      */
      Grid();

      /**
      * Default destructor.
      */
      virtual ~Grid();

      /**
      * Set Grid length and number of grid.
      */
      void setBox(const Vector& boxIn, const IntVector& nGridIn);

      /**
      * Rescale box dimension.
      *
      * \param s1  scaling factor for box edge length.
      */
      virtual void rescaleBox(const double s1);

      /**
      * Fill the list of site indices for a given vector.
      *
      * \param flag  > 0: particle insersion
      *              < 0: particle removal
      */
      void fillSiteList(const Vector& r, const int flag);

      /**
      * Find the list of sites affected by a inersting a particle at r.
      */
      void findAffectedSites(const Vector& r, const int dir = 2);

      /**
      * Find the list of sites affected by a virtual move and the corresponding weights.
      */
      void findAffectedSites(const Vector& r1, const Vector& r2);

      /**
      * Accessor to the number of changed sites.
      */
      int getNChanged() const;

      /**
      * Accessor to the list of changed site ids.
      */
      const std::vector< std::pair<int,double> >& getChangedSiteList() const;

      /**
      * Wrap 3d indices to 1d index.
      */
      int wrap(const IntVector& ind) const;

      /**
      * Unwrap 1d indices to 3d index.
      */
      void unwrap(const int id, IntVector& ind) const;

   protected:

      /// Box dimension.
      Vector boxL_;

      /// Box volume.
      double boxV_;

      /// Number of grid points.
      IntVector nGrid_;

      /// Half number of grid points along each direction.
      IntVector nGridHalf_;

      /// Grid size along each direction.
      Vector drGrid_;

      /// Lists of sites affected by a particle move.
      std::vector< std::pair<int,double> >   insertId_;
      std::vector< std::pair<int,double> >   removeId_;
      std::vector< std::pair<int,double> >   changeId_;
      int                                    nChanged_;

      /// Green's function for the Poisson's equation.
      PoissonGreenFunc poissonGreenFunc_;

      /**
      * Merge the insersion and removal lists.
      */
      void mergeLists();

      /**
      * Write the sum of sorted weights (diagnostic method).
      */
      void writeWeight(std::ostream& out);

      /**
      * Shift the index according to the minimum image convention.
      */
      void shiftDIndex(IntVector& dv) const;

   };

   /*
   * Wrap 3d indices to 1d index.
   */
   inline int Grid::wrap(const IntVector& ind) const
   { return ((ind[0]*nGrid_[1] + ind[1])*nGrid_[2] + ind[2]); }

   /*
   * Unwrap 1d indices to 3d index.
   */
   inline void Grid::unwrap(const int id, IntVector& ind) const
   {
      ind[2] = id % nGrid_[2];
      ind[1] = ((id - ind[2]) / nGrid_[2]) % nGrid_[1];
      ind[0] = (id - ind[2] - ind[1] * nGrid_[2]) / nGrid_[1] / nGrid_[2];
   }

   /**
   * Accessor to the number of changed sites.
   */
   inline int Grid::getNChanged() const
   { return nChanged_; }

   /**
   * Accessor to the list of changed site ids.
   */
   inline const std::vector< std::pair<int,double> >& Grid::getChangedSiteList() const
   { return changeId_; }

   /*
   * Shift the index difference according to the minimum image convention.
   */
   inline void Grid::shiftDIndex(IntVector& dv) const
   {
      for (int i = 0; i < Dimension; ++i) {
         if (dv[i] >= nGridHalf_[i])
            dv[i] -= nGrid_[i];
         else if (dv[i] <= -nGridHalf_[i])
            dv[i] += nGrid_[i];
      }
   }

}
#endif
