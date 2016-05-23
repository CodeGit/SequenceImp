/**
 * This file has no copyright assigned and is placed in the Public Domain.
 * This file is part of the w64 mingw-runtime package.
 * No warranty is given; refer to the file DISCLAIMER.PD within this package.
 */
#ifndef _AgentError_h_
#define _AgentError_h_

#define AGENTERROR(x) MAKE_SCODE(SEVERITY_ERROR,FACILITY_ITF,(x)+0x2000)
#define AGENTWARNING(x) MAKE_SCODE(SEVERITY_SUCCESS,FACILITY_ITF,(x)+0x2000)
#define AGENTREQERROR(x) MAKE_SCODE(SEVERITY_ERROR,FACILITY_ITF,(x)+0x2100)
#define AGENTPROVIDERERROR(x) MAKE_SCODE(SEVERITY_ERROR,FACILITY_ITF,(x)+0x2200)
#define AGENTVOICEERROR(x) MAKE_SCODE(SEVERITY_ERROR,FACILITY_ITF,(x)+0x2300)
#define AGENTAUDIOERROR(x) MAKE_SCODE(SEVERITY_ERROR,FACILITY_ITF,(x)+0x2400)
#define AGENTCTLERROR(x) MAKE_SCODE(SEVERITY_ERROR,FACILITY_ITF,(x)+0x2500)
#define AGENTEXTERROR(x) MAKE_SCODE(SEVERITY_ERROR,FACILITY_ITF,(x)+0x2600)

#define AGENTERR_CLIENTINVALID AGENTERROR(1)
#define AGENTERR_CHARACTERINVALID AGENTERROR(2)
#define AGENTERR_ANIMATIONNOTFOUND AGENTERROR(3)
#define AGENTERR_STATENOTFOUND AGENTERROR(4)
#define AGENTERR_AUDIONOTFOUND AGENTERROR(5)
#define AGENTERR_COMMANDNOTFOUND AGENTERROR(6)
#define AGENTERR_COMMANDALREADYINUSE AGENTERROR(7)
#define AGENTERR_MENUNOTFOUND AGENTERROR(8)
#define AGENTERR_LOSTCONNECTION AGENTERROR(9)
#define AGENTERR_CHARACTERNOTVISIBLE AGENTERROR(10)
#define AGENTERR_CHARACTERALREADYLOADED AGENTERROR(11)
#define AGENTERR_NOBALLOON AGENTERROR(12)
#define AGENTERR_NOCOMMANDSWINDOW AGENTERROR(13)
#define AGENTERR_INVALIDPREPARETYPE AGENTERROR(14)
#define AGENTERR_INVALIDANIMATION AGENTERROR(15)
#define AGENTERR_CANTMOVEDURINGDRAG AGENTERROR(16)
#define AGENTERR_CHARACTERNOTACTIVE AGENTERROR(17)
#define AGENTERR_LANGUAGENOTFOUND AGENTERROR(18)
#define AGENTERR_TTSLANGUAGENOTFOUND AGENTERROR(19)
#define AGENTERR_SRLANGUAGENOTFOUND AGENTERROR(20)
#define AGENTERR_LANGUAGEMISMATCH AGENTERROR(21)
#define AGENTERR_SPEAKINGDISABLED AGENTERROR(22)
#define AGENTERR_NOCHARACTERS AGENTERROR(23)
#define AGENTERR_DEFAULTCHARACTER AGENTERROR(24)

#define AGENTWARNING_TTSENGINENOTFOUND AGENTWARNING(1)
#define AGENTWARNING_ONLYCLIENT AGENTWARNING(2)

#define AGENTREQERR_OBJECTNOTFOUND AGENTREQERROR(1)
#define AGENTREQERR_OBJECTINVALID AGENTREQERROR(2)
#define AGENTREQERR_CANTSTOPOTHERS AGENTREQERROR(3)
#define AGENTREQERR_CANTINTERRUPTSELF AGENTREQERROR(4)
#define AGENTREQERR_CANTWAITONSELF AGENTREQERROR(5)
#define AGENTREQERR_INVALIDBOOKMARK AGENTREQERROR(6)
#define AGENTREQERR_SUSPENDED AGENTREQERROR(7)
#define AGENTREQERR_REMOVED AGENTREQERROR(8)

#define IS_INTERRUPT_ERROR(hRes) ((hRes >= AGENTREQERR_INTERRUPTEDLISTENKEY) && (hRes <= AGENTREQERR_INTERRUPTEDUSER))

#define AGENTREQERR_INTERRUPTEDLISTENKEY AGENTREQERROR(10)
#define AGENTREQERR_INTERRUPTEDHEARING AGENTREQERROR(11)
#define AGENTREQERR_INTERRUPTEDCODE AGENTREQERROR(12)
#define AGENTREQERR_INTERRUPTEDUSER AGENTREQERROR(13)

#define AGENTREQERR_INVALIDLASTTAG AGENTREQERROR(14)

#define AGENTPROVERROR_INIT AGENTPROVIDERERROR(1)
#define AGENTPROVERROR_CHARACTERVERSION AGENTPROVIDERERROR(2)
#define AGENTPROVERROR_VERSION AGENTPROVIDERERROR(3)
#define AGENTPROVERROR_MAGIC AGENTPROVIDERERROR(4)
#define AGENTPROVERROR_CHARACTERINVALID AGENTPROVIDERERROR(5)
#define AGENTPROVERROR_WAVEINVALID AGENTPROVIDERERROR(6)
#define AGENTPROVERROR_WAVECORRUPT AGENTPROVIDERERROR(7)
#define AGENTPROVERROR_MMIO AGENTPROVIDERERROR(8)
#define AGENTPROVERROR_PROTOCOL AGENTPROVIDERERROR(9)

#define AGENTAUDIOERROR_DEVICE AGENTAUDIOERROR(1)
#define AGENTAUDIOERROR_TTSENUMERATOR AGENTAUDIOERROR(2)
#define AGENTAUDIOERROR_TTSSELECT AGENTAUDIOERROR(3)
#define AGENTAUDIOERROR_TTSREGISTER AGENTAUDIOERROR(4)
#define AGENTAUDIOERROR_TTSUNEXPECTED AGENTAUDIOERROR(5)
#define AGENTAUDIOERROR_LWVINIT AGENTAUDIOERROR(6)
#define AGENTAUDIOERROR_LWVREGISTER AGENTAUDIOERROR(7)
#define AGENTAUDIOERROR_LWVUNEXPECTED AGENTAUDIOERROR(8)

#define AGENTCTLERROR_NOEVENTSAVAILABLE AGENTCTLERROR(1)
#define AGENTCTLERROR_SERVERINIT AGENTCTLERROR(2)
#define AGENTCTLERROR_LANGUAGE AGENTCTLERROR(3)

#define AGENTVOICEERROR_COULDNTSTARTDEVICE AGENTVOICEERROR(1)
#define AGENTVOICEERROR_NOTINSTALLED AGENTVOICEERROR(2)
#define AGENTVOICEERROR_NOTINITIALIZED AGENTVOICEERROR(3)
#define AGENTVOICEERROR_INVALIDMENU AGENTVOICEERROR(4)
#define AGENTVOICEERROR_UNCLOSEDALTERNATIVE AGENTVOICEERROR(5)
#define AGENTVOICEERROR_UNCLOSEDOPTIONAL AGENTVOICEERROR(6)
#define AGENTVOICEERROR_UNEXPECTEDENDOFALTERNATIVE AGENTVOICEERROR(7)
#define AGENTVOICEERROR_UNEXPECTEDENDOFOPTIONAL AGENTVOICEERROR(8)
#define AGENTVOICEERROR_UNEXPECTEDALTERNATIVE AGENTVOICEERROR(9)
#define AGENTVOICEERROR_NOSRMODE AGENTVOICEERROR(10)
#define AGENTVOICEERROR_SRMODENOTFOUND AGENTVOICEERROR(11)
#define AGENTVOICEERROR_SPEECHDISABLED AGENTVOICEERROR(12)
#define AGENTVOICEERROR_UNEXPECTEDENDOFREPEAT AGENTVOICEERROR(13)
#define AGENTVOICEERROR_UNCLOSEDREPEAT AGENTVOICEERROR(14)
#define AGENTVOICEERROR_UNEXPECTEDREPEAT AGENTVOICEERROR(15)

#define AGENTEXTERROR_EXTNOTFOUND AGENTEXTERROR(1)
#define AGENTEXTERROR_INVALIDCLIENT AGENTEXTERROR(2)
#endif
