#ifndef MPI_TRAITS_H
#define MPI_TRAITS_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>

#ifdef UTIL_MPI

namespace Util 
{

   /**
   * Default MpiTraits class.
   *
   * Each explicit specialization of MpiTraits has a public
   * static const member named type. This is the MPI data
   * type associated with the C++ template type parameter.
   */
   template <typename T>
   class MpiTraits
   { };

   /**
   * MpiTraits<char> explicit specialization.
   */
   template <>
   class MpiTraits<char>
   {  
   public: 
      static const MPI::Datatype type;   ///< MPI Datatype
      static const bool hasType;         ///< Is the MPI type initialized?
   };

   /**
   * MpiTraits<unsigned char> explicit specialization.
   */
   template <>
   class MpiTraits<unsigned char>
   {  
   public: 
      static const MPI::Datatype type;   ///< MPI Datatype
      static const bool hasType;         ///< Is the MPI type initialized?
   };

   /**
   * MpiTraits<short> explicit specialization.
   */
   template <>
   class MpiTraits<short>
   {  
   public: 
      static const MPI::Datatype type;   ///< MPI Datatype
      static const bool hasType;         ///< Is the MPI type initialized?
   };

   /**
   * MpiTraits<int> explicit specialization.
   */
   template <>
   class MpiTraits<int>
   {  
   public: 
      static const MPI::Datatype type;   ///< MPI Datatype
      static const bool hasType;         ///< Is the MPI type initialized?
   };

   /**
   * MpiTraits<long> explicit specialization.
   */
   template <>
   class MpiTraits<long>
   {  
   public: 
      static const MPI::Datatype type;   ///< MPI Datatype
      static const bool hasType;         ///< Is the MPI type initialized?
   };

   /**
   * MpiTraits<unsigned short> explicit specialization.
   */
   template <>
   class MpiTraits<unsigned short>
   {  
   public: 
      static const MPI::Datatype type;   ///< MPI Datatype
      static const bool hasType;         ///< Is the MPI type initialized?
   };

   /**
   * MpiTraits<unsigned int> explicit specialization.
   */
   template <>
   class MpiTraits<unsigned int>
   {  
   public: 
      static const MPI::Datatype type;  ///< MPI Datatype
      static const bool hasType;        ///< Is the MPI type initialized?
   };

   /**
   * MpiTraits<unsigned long> explicit specialization.
   */
   template <>
   class MpiTraits<unsigned long>
   {  
   public: 
      static const MPI::Datatype type;  ///< MPI Datatype
      static const bool hasType;        ///< Is the MPI type initialized?
   };

   /**
   * MpiTraits<float> explicit specialization.
   */
   template <>
   class MpiTraits<float>
   {  
   public: 
      static const MPI::Datatype type;  ///< MPI Datatype
      static const bool hasType;        ///< Is the MPI type initialized?
   };

   /**
   * MpiTraits<double> explicit specialization.
   */
   template <>
   class MpiTraits<double>
   {  
   public: 
      static const MPI::Datatype type; ///< MPI Datatype
      static const bool hasType;       ///< Is the MPI type initialized?
   };

   /**
   * MpiTraits<long double> explicit specialization.
   */
   template <>
   class MpiTraits<long double>
   {  
   public: 
      static const MPI::Datatype type;  ///< MPI Datatype 
      static const bool hasType;        ///< Is the MPI type initialized?
   };

   /**
   * MpiTraits<bool> explicit specialization.
   */
   template <>
   class MpiTraits<bool>
   {  
   public: 
      static const MPI::Datatype type;  ///< MPI Datatype
      static const bool hasType;        ///< Is the MPI type initialized?
   };

   #if 0
   /**
   * MpiTraits<wchar_t> explicit specialization.
   */
   template <>
   class MpiTraits<wchar_t>
   {  
   public: 
      static const MPI::Datatype type;  ///< MPI Datatype
      static const bool hasType;        ///< Is the MPI type initialized?
   };

   /**
   * MpiTraits< std::complex<float> > explicit specialization.
   */
   template <>
   class MpiTraits< std::complex<float> >
   {  
   public: 
      static const MPI::Datatype type;  ///< MPI Datatype
      static const bool hasType;        ///< Is the MPI type initialized?
   };

   /**
   * MpiTraits< std::complex<double> > explicit specialization.
   */
   template <>
   class MpiTraits< std::complex<double> >
   {  
   public: 
      static const MPI::Datatype type;  ///< MPI Datatype
      static const bool hasType;        ///< Is the MPI type initialized?
   };

   /**
   * MpiTraits< std::complex<long double> > explicit specialization.
   */
   template <>
   class MpiTraits< std::complex<long double> >
   {  
   public: 
      static const MPI::Datatype type; ///< MPI Datatype
      static const bool hasType;       ///< Is the MPI type initialized?
   };
   #endif

}
#endif
#endif
