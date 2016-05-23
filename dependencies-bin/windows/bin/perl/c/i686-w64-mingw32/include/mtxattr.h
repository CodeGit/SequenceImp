/**
 * This file has no copyright assigned and is placed in the Public Domain.
 * This file is part of the w64 mingw-runtime package.
 * No warranty is given; refer to the file DISCLAIMER.PD within this package.
 */
#ifndef _MTXATTR_H_
#define _MTXATTR_H_

#define TLBATTR_TRANS_REQUIRED 17093CC5-9BD2-11cf-AA4F-304BF89C0001
#define TLBATTR_TRANS_NOTSUPP 17093CC6-9BD2-11cf-AA4F-304BF89C0001
#define TLBATTR_TRANS_REQNEW 17093CC7-9BD2-11cf-AA4F-304BF89C0001
#define TLBATTR_TRANS_SUPPORTED 17093CC8-9BD2-11cf-AA4F-304BF89C0001
#define TLBATTR_QUEUEABLE E5FC3761-0BBA-11d2-B8FE-00C04FC340EE
#define TLBATTR_COMTI_INTRINSICS 47065EDC-D7FE-4B03-919C-C4A50B749605

#define TRANSACTION_REQUIRED custom(TLBATTR_TRANS_REQUIRED,0)
#define TRANSACTION_SUPPORTED custom(TLBATTR_TRANS_SUPPORTED,0)
#define TRANSACTION_NOT_SUPPORTED custom(TLBATTR_TRANS_NOTSUPP,0)
#define TRANSACTION_REQUIRES_NEW custom(TLBATTR_TRANS_REQNEW,0)
#define QUEUEABLE custom(TLBATTR_QUEUEABLE,0)
#define COMTI_INTRINSICS_ENABLED custom(TLBATTR_COMTI_INTRINSICS,0)

#endif
