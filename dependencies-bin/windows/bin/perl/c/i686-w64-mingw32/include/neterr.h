/**
 * This file has no copyright assigned and is placed in the Public Domain.
 * This file is part of the w64 mingw-runtime package.
 * No warranty is given; refer to the file DISCLAIMER.PD within this package.
 */
#ifndef _NMNETERR_
#define _NMNETERR_

#define NETERR_RING_STATUS_SIGNAL_LOST 0x00008000
#define NETERR_RING_STATUS_HARD_ERROR 0x00004000
#define NETERR_RING_STATUS_SOFT_ERROR 0x00002000
#define NETERR_RING_STATUS_TRANSMIT_BEACON 0x00001000
#define NETERR_RING_STATUS_LOBE_WIRE_FAULT 0x00000800
#define NETERR_RING_STATUS_AUTO_REMOVAL_ERROR 0x00000400
#define NETERR_RING_STATUS_REMOTE_RECEIVED 0x00000200
#define NETERR_RING_STATUS_COUNTER_OVERFLOW 0x00000100
#define NETERR_RING_STATUS_SIGNAL_STATION 0x00000080
#define NETERR_RING_STATUS_RECOVERY 0x00000040
#define NETERR_RING_STOP_CAPTURE 0x00008E00

#endif
