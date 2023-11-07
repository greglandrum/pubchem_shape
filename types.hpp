#ifndef ALIGN3D_TYPES__HPP
#define ALIGN3D_TYPES__HPP

#if (__SIZEOF_LONG__ == 8)
    typedef signed long Int8;    
    typedef unsigned long Uint8;
#elif (__SIZEOF_LONG_LONG__ == 8)
    typedef signed long long Int8;    
    typedef unsigned long long Uint8;
#else
    #error "This platform does not support 8-byte integer"
#endif

#endif // ALIGN3D_TYPES__HPP
