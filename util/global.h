#ifndef UTIL_GLOBAL_H
#define UTIL_GLOBAL_H

/**
* \file  util/global.h
* \brief File containing preprocessor macros for error handling.
*/

//-----------------------------------------------------------------------------
// MPI Parallelization

#ifdef UTIL_MPI
#include <mpi.h>
#endif

//-----------------------------------------------------------------------------
// Log file

/**
* Include access to a Log file.
*/
#include "Log.h"

//-----------------------------------------------------------------------------
// Errors: Assertions and and Exception macros

#ifndef UTIL_DEBUG

/**
* If defined, disable all C assert(...) statements. 
*/
#define NDEBUG        

#endif

#include <assert.h>
#include "Exception.h"

/**
* Macro for the name of the current function (compiler dependent).
*/
#define UTIL_FUNC __PRETTY_FUNCTION__

/**
* Macro for throwing an Exception, reporting function, file and line number.
*/
#ifdef  UTIL_FUNC
  #ifndef UTIL_MPI
    #define UTIL_THROW(msg) throw Exception(UTIL_FUNC, msg, __FILE__, __LINE__)
  #else
    #define UTIL_THROW(msg) { \
      Exception e(UTIL_FUNC, msg, __FILE__, __LINE__); \
      std::cerr   << e.message() << std::endl; \
      Log::file().flush(); \
      Log::close(); \
      MPI::COMM_WORLD.Abort(65); }
  #endif
#else
  #ifndef UTIL_MPI
    #define UTIL_THROW(msg) throw Exception(msg, __FILE__, __LINE__)
  #else
    #define UTIL_THROW(msg) { \
      Exception e(msg, __FILE__, __LINE__); \
      std::cerr   << e.message() << std::endl; \
      Log::file().flush(); \
      Log::close(); \
      MPI::COMM_WORLD.Abort(65); }
  #endif
#endif

#endif
