#pragma once

#ifdef AADC_DLL
    #ifdef AADC_EXPORTS
        #define AADC_API __declspec(dllexport)
    #else
        #define AADC_API __declspec(dllimport)
    #endif // ifdef AADC_EXPORTS
#else
    #define AADC_API
#endif // ifdef AADC_DLL