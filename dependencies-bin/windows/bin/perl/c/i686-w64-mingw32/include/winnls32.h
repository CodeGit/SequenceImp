/**
 * This file has no copyright assigned and is placed in the Public Domain.
 * This file is part of the w64 mingw-runtime package.
 * No warranty is given; refer to the file DISCLAIMER.PD within this package.
 */
#ifndef _WINNLS32_
#define _WINNLS32_

#include <_mingw_unicode.h>

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct _tagDATETIME {
    WORD year;
    WORD month;
    WORD day;
    WORD hour;
    WORD min;
    WORD sec;
  } DATETIME;

  typedef struct _tagIMEPROA {
    HWND hWnd;
    DATETIME InstDate;
    UINT wVersion;
    BYTE szDescription[50];
    BYTE szName[80];
    BYTE szOptions[30];
  } IMEPROA,*PIMEPROA,NEAR *NPIMEPROA,*LPIMEPROA;

  typedef struct _tagIMEPROW {
    HWND hWnd;
    DATETIME InstDate;
    UINT wVersion;
    WCHAR szDescription[50];
    WCHAR szName[80];
    WCHAR szOptions[30];
  } IMEPROW,*PIMEPROW,NEAR *NPIMEPROW,*LPIMEPROW;

  __MINGW_TYPEDEF_AW(IMEPRO)
  __MINGW_TYPEDEF_AW(PIMEPRO)
  __MINGW_TYPEDEF_AW(NPIMEPRO)
  __MINGW_TYPEDEF_AW(LPIMEPRO)

#define IMPGetIME __MINGW_NAME_AW(IMPGetIME)
#define IMPQueryIME __MINGW_NAME_AW(IMPQueryIME)
#define IMPSetIME __MINGW_NAME_AW(IMPSetIME)

  WINBOOL WINAPI IMPGetIMEA(HWND,LPIMEPROA);
  WINBOOL WINAPI IMPGetIMEW(HWND,LPIMEPROW);
  WINBOOL WINAPI IMPQueryIMEA(LPIMEPROA);
  WINBOOL WINAPI IMPQueryIMEW(LPIMEPROW);
  WINBOOL WINAPI IMPSetIMEA(HWND,LPIMEPROA);
  WINBOOL WINAPI IMPSetIMEW(HWND,LPIMEPROW);
  UINT WINAPI WINNLSGetIMEHotkey(HWND);
  WINBOOL WINAPI WINNLSEnableIME(HWND,WINBOOL);
  WINBOOL WINAPI WINNLSGetEnableStatus(HWND);

#ifdef __cplusplus
}
#endif
#endif
