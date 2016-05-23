/**
 * This file has no copyright assigned and is placed in the Public Domain.
 * This file is part of the w64 mingw-runtime package.
 * No warranty is given; refer to the file DISCLAIMER.PD within this package.
 */

#ifndef __WPCEVENT_H__
#define __WPCEVENT_H__

#include <evntprov.h>

EXTERN_C DECLSPEC_SELECTANY const GUID WPCPROV = {0x01090065, 0xb467, 0x4503, {0x9b,0x28,0x53,0x37,0x66,0x76,0x10,0x87}};

typedef enum tagWPC_ARGS_FILEDOWNLOADEVENT {
    WPC_ARGS_FILEDOWNLOADEVENT_URL = 0,
    WPC_ARGS_FILEDOWNLOADEVENT_APPNAME,
    WPC_ARGS_FILEDOWNLOADEVENT_VERSION,
    WPC_ARGS_FILEDOWNLOADEVENT_BLOCKED,
    WPC_ARGS_FILEDOWNLOADEVENT_PATH,
    WPC_ARGS_FILEDOWNLOADEVENT_CARGS
} WPC_ARGS_FILEDOWNLOADEVENT;

EXTERN_C DECLSPEC_SELECTANY const EVENT_DESCRIPTOR WPCEVENT_WEB_FILEDOWNLOAD = {0xa,0x0,0x10,0x4,0x18,0xa,0x8000000000000030};

#define WPCEVENT_WEB_FILEDOWNLOAD_value 0xa

#endif /*__WPCEVENT_H__*/
