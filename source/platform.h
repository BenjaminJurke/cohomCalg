////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  platform.h                                                       +------------------------+   //
//  ==========                                                       | Cross platform handler |   //
//                                                                   +------------------------+   //
//  Code: Benjamin Jurke, http://benjaminjurke.net                                                //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//  File history:                                                                                 //
//        - 27.04.2010  File created as platform.h                                                //
//                      This file takes care of a couple of definitions and environment variables //
//                      which are compiler and/or platform specific.                              //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef INC_PLATFORM_H
#define INC_PLATFORM_H


// NOTE: I've only checked this code on MSVS2010 and Ubuntu/gcc.


// We define the command save_sprintf which is the ordinary sprintf string formatting command
// but with an additional buffer length check which avoids the buffer overflow errors / security leaks
// quite a number of C/C++ programs suffer from.
//
// Following the C99 standard, a new function snprintf was introduced to take care of this issue.
// However, for god-knows-why Microsoft decided to use a different name, so we are using a define to
// allow for further ingenious cross-platform problems with this...
//
// NOTE: You still have to include <cstring> to use this function...

#if defined(_MSC_VER)                   // MS compiler (Windows)
    #define safe_sprintf sprintf_s
    #define string_to_int64(__v_val_) _atoi64(__v_val_)

	#include <Windows.h>
	inline void SleepMilliSec(unsigned long ulMilliseconds)
	{
		Sleep(ulMilliseconds);
	}

#elif defined (__GNUC__)                // gcc compiler (Unix/Linux/MacOS)

    #define safe_sprintf snprintf
    #define string_to_int64(__v_val_) atoll(__v_val_)
    #include <time.h>
    inline void SleepMilliSec(unsigned long ulMilliseconds)
    {
        struct timespec req = {0};
        req.tv_sec = 0;
        req.tv_nsec = ulMilliseconds * 1000000L;
        nanosleep(&req, (struct timespec *) NULL);
    }

#else
    #error Unrecognized compiler.
#endif


// Next we have a couple of platform variables, which are just for convenience output

#if defined(_WIN32)
    #define APP_TARGET_OS "Windows"
#elif defined(__unix__)
    #define APP_TARGET_OS "Linux/Unix"
#elif (defined(__APPLE__) && defined(__MACH__))
    #define APP_TARGET_OS "MacOS X"
#else
    #define APP_TARGET_OS "Unrecognized OS"
#endif

#if (defined(_M_X64) || defined(__amd64__))
    #define APP_ARCH "x86-64 / 64 bit"
#elif (defined(_M_IA64))
    #define APP_ARCH "ia64 (Itanium) / 64 bit"
#elif (defined(_M_IX86) || defined (__i386__))
    #define APP_ARCH "i386 / 32 bit"
#else
    #define APP_ARCH "Unrecognized"
#endif


#endif
