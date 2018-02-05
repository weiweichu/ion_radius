#ifndef CONSTANTS_H
#define CONSTANTS_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <complex>

namespace Util
{

   /**
   * Mathematical constants.
   *
   * \ingroup Math_Module
   */
   class Constants
   {

   public:

      /**
      * Initialize static constants.
      */
      static void initStatic();

      /**
      * Square root of -1.
      */
      static const double               Epsilon;

      /**
      * Trigonometric constant Pi.
      */
      static const double               Pi;

      /**
      * Square root of -1.
      */
      static const std::complex<double> Im;

   };

} 
#endif
