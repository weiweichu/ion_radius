#ifndef GRID_CPP
#define GRID_CPP

/*
* Copyright 2010, Jian Qin (qinjian000@gmail.com)
* Distributed under the terms of the GNU General Public License.
*/

#include "Grid.h"
#include "../util/Constants.h"

#include <math.h>

namespace GridMC
{
   using namespace std;

   /*
   * Constructor.
   *
   * Assumptions:
   *    (1) Unit cell shape is orthogonal.
   *    (2) PM1 scheme is used, so the dimensions 8, 8, and 16.
   */
   Grid::Grid():
      boxL_(1.0),
      boxV_(1.0),
      nGrid_(2),
      nGridHalf_(1),
      drGrid_(0.5),
      poissonGreenFunc_()
   {
      std::pair<int,double> zero(0,0.0);
      insertId_.resize(8, zero);
      removeId_.resize(8, zero);
      changeId_.resize(16, zero);
      nChanged_ = 0;
   }

   /*
   * Default destructor.
   */
   Grid::~Grid()
   {}

   /*
   * Set unit cell length and grid dimension.
   */
   void Grid::setBox(const Vector& boxLIn, const IntVector& nGridIn)
   {
      int i;

      boxL_ = boxLIn;
      boxV_ = 1.0;
      for (i = 0; i < Dimension; ++i) boxV_ *= boxL_[i];

      nGrid_ = nGridIn;
      for (i = 0; i < Dimension; ++i) {
         if (nGrid_[i]%2 == 0)
            nGridHalf_[i] = nGrid_[i] / 2 + 1;
         else
            nGridHalf_[i] = (nGrid_[i] + 1) / 2;
      }

      for (i = 0; i < Dimension; ++i)
         drGrid_[i] = boxL_[i] / double(nGrid_[i]);

      poissonGreenFunc_.setBox(boxLIn, nGridIn);
   }


   /*
   * Reset unit cell length and grid dimension.
   */
   void Grid::rescaleBox(const double s1)
   {
      boxL_   *= s1;
      boxV_   *= s1 * s1 * s1;
      drGrid_ *= s1;

      poissonGreenFunc_.rescaleBox(s1);
   }


   /*
   * Fill the list of indices affected a particle move.
   *
   * Indices are ordered by the order of increasing id.
   */
   void Grid::fillSiteList(const Vector& r, const int flag)
   {
      int    i, ic, ix, iy, iz;
      int    idx[2], idy[2], idz[3];
      double wtx[2], wty[2], wtz[3];
      IntVector ind0, ind;
      Vector rCell;

      // Find the primary cell indices and relative weight.
      for (i = 0; i < 3; ++i) {
         rCell[i] = r[i] - floor(r[i] / boxL_[i]) * boxL_[i];
         ind0[i] = int(floor(rCell[i] / drGrid_[i]));
         rCell[i] -= drGrid_[i] * double(ind0[i]);
      }

      // Order the cell by id.
      if (ind0[0] < nGrid_[0] - 1) {
         idx[0] = ind0[0];
         idx[1] = ind0[0] + 1;
         wtx[1] = rCell[0] / drGrid_[0];
         wtx[0] = 1.0 - wtx[1];
      } else {
         idx[1] = ind0[0];
         idx[0] = 0;
         wtx[0] = rCell[0] / drGrid_[0];
         wtx[1] = 1.0 - wtx[0];
      }

      if (ind0[1] < nGrid_[1] - 1) {
         idy[0] = ind0[1];
         idy[1] = ind0[1] + 1;
         wty[1] = rCell[1] / drGrid_[1];
         wty[0] = 1.0 - wty[1];
      } else {
         idy[1] = ind0[1];
         idy[0] = 0;
         wty[0] = rCell[1] / drGrid_[1];
         wty[1] = 1.0 - wty[0];
      }

      if (ind0[2] < nGrid_[2] - 1) {
         idz[0] = ind0[2];
         idz[1] = ind0[2] + 1;
         wtz[1] = rCell[2] / drGrid_[2];
         wtz[0] = 1.0 - wtz[1];
      } else {
         idz[1] = ind0[2];
         idz[0] = 0;
         wtz[0] = rCell[2] / drGrid_[2];
         wtz[1] = 1.0 - wtz[0];
      }

      // Fill index array.
      ic = 0;
      if (flag == 1) {

         for (ix = 0; ix < 2; ++ix) {
         for (iy = 0; iy < 2; ++iy) {
         for (iz = 0; iz < 2; ++iz) {
            ind[0] = idx[ix];
            ind[1] = idy[iy];
            ind[2] = idz[iz];
            insertId_[ic].first  = wrap(ind);
            insertId_[ic++].second = wtx[ix] * wty[iy] * wtz[iz];
         } } }

      } else if (flag == -1) {

         for (ix = 0; ix < 2; ++ix) {
         for (iy = 0; iy < 2; ++iy) {
         for (iz = 0; iz < 2; ++iz) {
            ind[0] = idx[ix];
            ind[1] = idy[iy];
            ind[2] = idz[iz];
            removeId_[ic].first  = wrap(ind);
            removeId_[ic++].second = -wtx[ix] * wty[iy] * wtz[iz];
         } } }

      } else if (flag == 2) {

         nChanged_ = 0;
         for (ix = 0; ix < 2; ++ix) {
         for (iy = 0; iy < 2; ++iy) {
         for (iz = 0; iz < 2; ++iz) {
            ind[0] = idx[ix];
            ind[1] = idy[iy];
            ind[2] = idz[iz];
            changeId_[nChanged_].first  = wrap(ind);
            changeId_[nChanged_++].second = wtx[ix] * wty[iy] * wtz[iz];
         } } }

      } else if (flag == -2) {

         nChanged_ = 0;
         for (ix = 0; ix < 2; ++ix) {
         for (iy = 0; iy < 2; ++iy) {
         for (iz = 0; iz < 2; ++iz) {
            ind[0] = idx[ix];
            ind[1] = idy[iy];
            ind[2] = idz[iz];
            changeId_[nChanged_].first  = wrap(ind);
            changeId_[nChanged_++].second = -wtx[ix] * wty[iy] * wtz[iz];
         } } }

      } else {

         UTIL_THROW("Invalid flag value.");

      }

   }

   /*
   * Merge the insertId and removeId lists.
   *
   * Assumption: insertId and removeId are ascending.
   */
   void Grid::mergeLists()
   {
      int i(0), j(0);

      nChanged_ = 0;
      while (i < 8 && j < 8) {
         if (insertId_[i].first < removeId_[j].first) {
            changeId_[nChanged_++] = insertId_[i++];
         } else if (insertId_[i].first > removeId_[j].first) {
            changeId_[nChanged_++] = removeId_[j++];
         } else {
            changeId_[nChanged_] = insertId_[i++];
            changeId_[nChanged_++].second += removeId_[j++].second;
         }
      }

      if (i == 8 && j < 8) {
         for (; j < 8; ++j) changeId_[nChanged_++] = removeId_[j];
      } else if (j == 8 && i < 8) {
         for (; i < 8; ++i) changeId_[nChanged_++] = insertId_[i];
      }
   }

   /*
   * Find the list of sites affected by a placing a particle at r.
   */
   void Grid::findAffectedSites(const Vector& r, const int dir)
   {
      fillSiteList(r, dir);
   }

   /*
   * Find the list of sites affected by a proposed r1 --> r2 move.
   */
   void Grid::findAffectedSites(const Vector& r1, const Vector& r2)
   {
      fillSiteList(r1, -1);
      fillSiteList(r2, +1);
      mergeLists();
   }

   /*
   * Output the summed weights for sorted arrays.
   */
   void Grid::writeWeight(std::ostream& out)
   {
      double w1(0.0), w2(0.0), w3(0.0);
      int i;

      out << "insertId: ";
      for (i = 0; i < 8; ++i) {
         out << "  " << insertId_[i].first;
         w1 += insertId_[i].second;
      }
      out << endl;

      out << "removeId: ";
      for (i = 0; i < 8; ++i) {
         out << "  " << removeId_[i].first;
         w2 += insertId_[i].second;
      }
      out << endl;

      out << "changeId: ";
      for (i = 0; i < nChanged_; ++i) {
         out << "  " << changeId_[i].first;
         w3 += changeId_[i].second;
      }
      out << endl;

      out << "sum(insert weight): " << w1;
      out << "   sum(remove weight): " << w2;
      out << "   sum(change weight): " << w3 << endl;
   }

}

#endif
