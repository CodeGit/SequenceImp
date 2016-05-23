/**
 * This file has no copyright assigned and is placed in the Public Domain.
 * This file is part of the w64 mingw-runtime package.
 * No warranty is given; refer to the file DISCLAIMER.PD within this package.
 */
#ifndef _DHCPCSDK_
#define _DHCPCSDK_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef DHCP_OPTIONS_DEFINED
#define DHCP_OPTIONS_DEFINED

#define OPTION_PAD 0
#define OPTION_SUBNET_MASK 1
#define OPTION_TIME_OFFSET 2
#define OPTION_ROUTER_ADDRESS 3
#define OPTION_TIME_SERVERS 4
#define OPTION_IEN116_NAME_SERVERS 5
#define OPTION_DOMAIN_NAME_SERVERS 6
#define OPTION_LOG_SERVERS 7
#define OPTION_COOKIE_SERVERS 8
#define OPTION_LPR_SERVERS 9
#define OPTION_IMPRESS_SERVERS 10
#define OPTION_RLP_SERVERS 11
#define OPTION_HOST_NAME 12
#define OPTION_BOOT_FILE_SIZE 13
#define OPTION_MERIT_DUMP_FILE 14
#define OPTION_DOMAIN_NAME 15
#define OPTION_SWAP_SERVER 16
#define OPTION_ROOT_DISK 17
#define OPTION_EXTENSIONS_PATH 18

#define OPTION_BE_A_ROUTER 19
#define OPTION_NON_LOCAL_SOURCE_ROUTING 20
#define OPTION_POLICY_FILTER_FOR_NLSR 21
#define OPTION_MAX_REASSEMBLY_SIZE 22
#define OPTION_DEFAULT_TTL 23
#define OPTION_PMTU_AGING_TIMEOUT 24
#define OPTION_PMTU_PLATEAU_TABLE 25

#define OPTION_MTU 26
#define OPTION_ALL_SUBNETS_MTU 27
#define OPTION_BROADCAST_ADDRESS 28
#define OPTION_PERFORM_MASK_DISCOVERY 29
#define OPTION_BE_A_MASK_SUPPLIER 30
#define OPTION_PERFORM_ROUTER_DISCOVERY 31
#define OPTION_ROUTER_SOLICITATION_ADDR 32
#define OPTION_STATIC_ROUTES 33
#define OPTION_TRAILERS 34
#define OPTION_ARP_CACHE_TIMEOUT 35
#define OPTION_ETHERNET_ENCAPSULATION 36

#define OPTION_TTL 37
#define OPTION_KEEP_ALIVE_INTERVAL 38
#define OPTION_KEEP_ALIVE_DATA_SIZE 39

#define OPTION_NETWORK_INFO_SERVICE_DOM 40
#define OPTION_NETWORK_INFO_SERVERS 41
#define OPTION_NETWORK_TIME_SERVERS 42

#define OPTION_VENDOR_SPEC_INFO 43

#define OPTION_NETBIOS_NAME_SERVER 44
#define OPTION_NETBIOS_DATAGRAM_SERVER 45
#define OPTION_NETBIOS_NODE_TYPE 46
#define OPTION_NETBIOS_SCOPE_OPTION 47

#define OPTION_XWINDOW_FONT_SERVER 48
#define OPTION_XWINDOW_DISPLAY_MANAGER 49

#define OPTION_REQUESTED_ADDRESS 50
#define OPTION_LEASE_TIME 51
#define OPTION_OK_TO_OVERLAY 52
#define OPTION_MESSAGE_TYPE 53
#define OPTION_SERVER_IDENTIFIER 54
#define OPTION_PARAMETER_REQUEST_LIST 55
#define OPTION_MESSAGE 56
#define OPTION_MESSAGE_LENGTH 57
#define OPTION_RENEWAL_TIME 58
#define OPTION_REBIND_TIME 59
#define OPTION_CLIENT_CLASS_INFO 60
#define OPTION_CLIENT_ID 61

#define OPTION_TFTP_SERVER_NAME 66
#define OPTION_BOOTFILE_NAME 67

#define OPTION_END 255
#endif

#ifndef DHCPAPI_PARAMS_DEFINED
#define DHCPAPI_PARAMS_DEFINED
  typedef struct _DHCPAPI_PARAMS {
    ULONG Flags;
    ULONG OptionId;
    WINBOOL IsVendor;
    LPBYTE Data;
    DWORD nBytesData;
  } DHCPAPI_PARAMS,*PDHCPAPI_PARAMS,*LPDHCPAPI_PARAMS;
#endif

  typedef struct _DHCPAPI_PARAMS DHCPCAPI_PARAMS,*PDHCPCAPI_PARAMS,*LPDHCPCAPI_PARAMS;

  typedef struct _DHCPCAPI_PARAMS_ARARAY {
    ULONG nParams;
    LPDHCPCAPI_PARAMS Params;
  } DHCPCAPI_PARAMS_ARRAY,*PDHCPCAPI_PARAMS_ARRAY,*LPDHCPCAPI_PARAMS_ARRAY;

  typedef struct _DHCPCAPI_CLASSID {
    ULONG Flags;
    LPBYTE Data;
    ULONG nBytesData;
  } DHCPCAPI_CLASSID,*PDHCPCAPI_CLASSID,*LPDHCPCAPI_CLASSID;

#define DHCPCAPI_REQUEST_PERSISTENT 0x01
#define DHCPCAPI_REQUEST_SYNCHRONOUS 0x02
#define DHCPCAPI_REQUEST_ASYNCHRONOUS 0x04
#define DHCPCAPI_REQUEST_CANCEL 0x08
#define DHCPCAPI_REQUEST_MASK 0x0F

  DWORD WINAPI DhcpCApiInitialize(LPDWORD Version);
  VOID WINAPI DhcpCApiCleanup(VOID);
  DWORD WINAPI DhcpRequestParams(DWORD Flags,LPVOID Reserved,LPWSTR AdapterName,LPDHCPCAPI_CLASSID ClassId,DHCPCAPI_PARAMS_ARRAY SendParams,DHCPCAPI_PARAMS_ARRAY RecdParams,LPBYTE Buffer,LPDWORD pSize,LPWSTR RequestIdStr);
  DWORD WINAPI DhcpUndoRequestParams(DWORD Flags,LPVOID Reserved,LPWSTR AdapterName,LPWSTR RequestIdStr);

#define DHCPCAPI_REGISTER_HANDLE_EVENT 0x01
  DWORD WINAPI DhcpRegisterParamChange(DWORD Flags,LPVOID Reserved,LPWSTR AdapterName,LPDHCPCAPI_CLASSID ClassId,DHCPCAPI_PARAMS_ARRAY Params,LPVOID Handle);

#define DHCPCAPI_DEREGISTER_HANDLE_EVENT 0x01
  DWORD WINAPI DhcpDeRegisterParamChange(DWORD Flags,LPVOID Reserved,LPVOID Event);
  DWORD WINAPI DhcpRemoveDNSRegistrations(VOID);

#ifdef __cplusplus
}
#endif
#endif
