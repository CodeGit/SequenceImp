/**
 * This file has no copyright assigned and is placed in the Public Domain.
 * This file is part of the w64 mingw-runtime package.
 * No warranty is given; refer to the file DISCLAIMER.PD within this package.
 */
#define categoryHeader 0x00000001L
#define categoryUnused 0x00000002L
#define categoryGeneral 0x00000003L
#define CDO_E_UNCAUGHT_EXCEPTION 0x80040201L
#define CDO_E_NOT_OPENED 0x80040202L
#define CDO_E_UNSUPPORTED_DATASOURCE 0x80040203L
#define CDO_E_INVALID_PROPERTYNAME 0x80040204L
#define CDO_E_PROP_UNSUPPORTED 0x80040205L
#define CDO_E_INACTIVE 0x80040206L
#define CDO_E_NO_SUPPORT_FOR_OBJECTS 0x80040207L
#define CDO_E_NOT_AVAILABLE 0x80040208L
#define CDO_E_NO_DEFAULT_DROP_DIR 0x80040209L
#define CDO_E_SMTP_SERVER_REQUIRED 0x8004020AL
#define CDO_E_NNTP_SERVER_REQUIRED 0x8004020BL
#define CDO_E_RECIPIENT_MISSING 0x8004020CL
#define CDO_E_FROM_MISSING 0x8004020DL
#define CDO_E_SENDER_REJECTED 0x8004020EL
#define CDO_E_RECIPIENTS_REJECTED 0x8004020FL
#define CDO_E_NNTP_POST_FAILED 0x80040210L
#define CDO_E_SMTP_SEND_FAILED 0x80040211L
#define CDO_E_CONNECTION_DROPPED 0x80040212L
#define CDO_E_FAILED_TO_CONNECT 0x80040213L
#define CDO_E_INVALID_POST 0x80040214L
#define CDO_E_AUTHENTICATION_FAILURE 0x80040215L
#define CDO_E_INVALID_CONTENT_TYPE 0x80040216L
#define CDO_E_LOGON_FAILURE 0x80040217L
#define CDO_E_HTTP_NOT_FOUND 0x80040218L
#define CDO_E_HTTP_FORBIDDEN 0x80040219L
#define CDO_E_HTTP_FAILED 0x8004021AL
#define CDO_E_MULTIPART_NO_DATA 0x8004021BL
#define CDO_E_INVALID_ENCODING_FOR_MULTIPART 0x8004021CL
#define CDO_E_UNSAFE_OPERATION 0x8004021DL
#define CDO_E_PROP_NOT_FOUND 0x8004021EL
#define CDO_E_INVALID_SEND_OPTION 0x80040220L
#define CDO_E_INVALID_POST_OPTION 0x80040221L
#define CDO_E_NO_PICKUP_DIR 0x80040222L
#define CDO_E_NOT_ALL_DELETED 0x80040223L
#define CDO_E_NO_METHOD 0x80040224L
#define CDO_E_PROP_READONLY 0x80040227L
#define CDO_E_PROP_CANNOT_DELETE 0x80040228L
#define CDO_E_BAD_DATA 0x80040229L
#define CDO_E_PROP_NONHEADER 0x8004022AL
#define CDO_E_INVALID_CHARSET 0x8004022BL
#define CDO_E_ADOSTREAM_NOT_BOUND 0x8004022CL
#define CDO_E_CONTENTPROPXML_NOT_FOUND 0x8004022DL
#define CDO_E_CONTENTPROPXML_WRONG_CHARSET 0x8004022EL
#define CDO_E_CONTENTPROPXML_PARSE_FAILED 0x8004022FL
#define CDO_E_CONTENTPROPXML_CONVERT_FAILED 0x80040230L
#define CDO_E_NO_DIRECTORIES_SPECIFIED 0x80040231L
#define CDO_E_DIRECTORIES_UNREACHABLE 0x80040232L
#define CDO_E_BAD_SENDER 0x80040233L
#define CDO_E_SELF_BINDING 0x80040234L
#define CDO_E_BAD_ATTENDEE_DATA 0x80040235L
#define CDO_E_ARGUMENT1 0x80044000L
#define CDO_E_ARGUMENT2 0x80044001L
#define CDO_E_ARGUMENT3 0x80044002L
#define CDO_E_ARGUMENT4 0x80044003L
#define CDO_E_ARGUMENT5 0x80044004L
#define CDO_E_NOT_FOUND 0x800CCE05L
#define CDO_E_INVALID_ENCODING_TYPE 0x800CCE1DL
#define IDS_ORIGINAL_MESSAGE 0x00011000L
#define IDS_FROM 0x00011001L
#define IDS_SENT 0x00011002L
#define IDS_POSTED_AT 0x00011003L
#define IDS_TO 0x00011004L
#define IDS_CC 0x00011005L
#define IDS_POSTED_TO 0x00011006L
#define IDS_CONVERSATION 0x00011007L
#define IDS_SUBJECT 0x00011008L
#define IDS_IMPORTANCE 0x00011009L
#define IDS_ON_BEHALF_OF 0x0001100AL
#define IDS_FW 0x0001100BL
#define IDS_RE 0x0001100CL
#define IDS_CODEPAGE 0x0001100DL

#ifdef CDOSVR
#define IDS_CalendarFolder 0x0001100EL
#define IDS_ContactsFolder 0x0001100FL
#define IDS_DraftsFolder 0x00011010L
#define IDS_JournalFolder 0x00011011L
#define IDS_NotesFolder 0x00011012L
#define IDS_TasksFolder 0x00011013L
#endif

#define IDS_NewFolder 0x00011014L
#define IDS_Location 0x00011015L
#define IDS_StartTime 0x00011016L
#define IDS_EndTime 0x00011017L
#define IDS_TimeZone 0x00011018L
#define IDS_LocalTime 0x00011019L
#define IDS_Organizer 0x0001101AL
#define IDS_ApptType 0x0001101BL
#define IDS_SingleAppt 0x0001101CL
#define IDS_SingleMtg 0x0001101DL
#define IDS_RecurAppt 0x0001101EL
#define IDS_RecurMtg 0x0001101FL
#define IDS_Universal 0x00011100L
#define IDS_Greenwich 0x00011101L
#define IDS_Sarajevo 0x00011102L
#define IDS_Paris 0x00011103L
#define IDS_Berlin 0x00011104L
#define IDS_EasternEurope 0x00011105L
#define IDS_Prague 0x00011106L
#define IDS_Athens 0x00011107L
#define IDS_Brasilia 0x00011108L
#define IDS_Atlantic 0x00011109L
#define IDS_Eastern 0x0001110AL
#define IDS_Central 0x0001110BL
#define IDS_Mountain 0x0001110CL
#define IDS_Pacific 0x0001110DL
#define IDS_Alaska 0x0001110EL
#define IDS_Hawaii 0x0001110FL
#define IDS_Midway 0x00011110L
#define IDS_Wellington 0x00011111L
#define IDS_Brisbane 0x00011112L
#define IDS_Adelaide 0x00011113L
#define IDS_Tokyo 0x00011114L
#define IDS_Singapore 0x00011115L
#define IDS_Bangkok 0x00011116L
#define IDS_Bombay 0x00011117L
#define IDS_AbuDhabi 0x00011118L
#define IDS_Tehran 0x00011119L
#define IDS_Baghdad 0x0001111AL
#define IDS_Israel 0x0001111BL
#define IDS_Newfoundland 0x0001111CL
#define IDS_Azores 0x0001111DL
#define IDS_MidAtlantic 0x0001111EL
#define IDS_Monrovia 0x0001111FL
#define IDS_BuenosAires 0x00011120L
#define IDS_Caracas 0x00011121L
#define IDS_Indiana 0x00011122L
#define IDS_Bogota 0x00011123L
#define IDS_Saskatchewan 0x00011124L
#define IDS_Mexico 0x00011125L
#define IDS_Arizona 0x00011126L
#define IDS_Eniwetok 0x00011127L
#define IDS_Fiji 0x00011128L
#define IDS_Magadan 0x00011129L
#define IDS_Hobart 0x0001112AL
#define IDS_Guam 0x0001112BL
#define IDS_Darwin 0x0001112CL
#define IDS_Beijing 0x0001112DL
#define IDS_Almaty 0x0001112EL
#define IDS_Islamabad 0x0001112FL
#define IDS_Kabul 0x00011130L
#define IDS_Cairo 0x00011131L
#define IDS_Harare 0x00011132L
#define IDS_Moscow 0x00011133L
#define IDS_CapeVerde 0x00011134L
#define IDS_Caucasus 0x00011135L
#define IDS_CentralAmerica 0x00011136L
#define IDS_EastAfrica 0x00011137L
#define IDS_Melbourne 0x00011138L
#define IDS_Ekaterinburg 0x00011139L
#define IDS_Helsinki 0x0001113AL
#define IDS_Greenland 0x0001113BL
#define IDS_Rangoon 0x0001113CL
#define IDS_Nepal 0x0001113DL
#define IDS_Irkutsk 0x0001113EL
#define IDS_Krasnoyarsk 0x0001113FL
#define IDS_Santiago 0x00011140L
#define IDS_SriLanka 0x00011141L
#define IDS_Tonga 0x00011142L
#define IDS_Vladivostok 0x00011143L
#define IDS_WestCentralAfrica 0x00011144L
#define IDS_Yakutsk 0x00011145L
#define IDS_Dhaka 0x00011146L
#define IDS_Seoul 0x00011147L
#define IDS_Perth 0x00011148L
#define IDS_Arab 0x00011149L
#define IDS_Taipei 0x0001114AL
#define IDS_Sydney2000 0x0001114BL

#ifdef CDOSVR
#define evtMethodCalled 0x00032000L
#define evtMethodReturning 0x00032001L
#define evtIsAborting 0xC0032002L
#define evtExpansionInitialized 0x00032003L
#define evtExpansionUnInitialized 0x00032004L
#define evtExpansionInitializeFailed 0xC0032005L
#define evtExpansionRegisterFailed 0xC0032006L
#define evtExpansionMessageSaveChangesFailed 0xC0032007L
#define evtExpansionMessageDeleteFailed 0xC0032008L
#define evtExpansionFolderSaveChangesFailed 0xC0032009L
#define evtExpansionTooManyInstancesPerDay 0x8003200AL
#define evtMailboxCreateTotalFailure 0xC003200BL
#define evtMailboxCreatePartialFailure 0xC003200CL
#define evtUninitImplRestFailed 0xC003200DL
#define evtExpandSavingAppt 0xC003200EL
#define evtExpandDeletingAppt 0xC003200FL
#define evtExpandQuery 0xC0032010L
#define evtExpandFolderSetProps 0xC0032011L
#define evtRegistryFailure 0xC0032012L
#define evtExpStat 0xC0032013L
#define evtDumpFcn 0xC0032014L
#define evtSaveDeleteFailFBUpdate 0xC0032015L
#define evtProcessingQueryCallback 0xC0032016L
#define evtMailboxLocalizeTotalFailure 0xC0032017L
#define evtMailboxLocalizePartialFailure 0xC0032018L
#define evtExpandMaster 0xC0032019L
#define evtExpansionInit 0xC003201AL
#define evtFBGenerateMsg 0xC003201BL
#define evtExpansionInstExpiryInPublicMDB 0x8003201CL
#define evtUnhandledExceptionInitialization 0xC003201DL
#define evtUnhandledExceptionShutdown 0xC003201EL
#define evtUnhandledExceptionInitializationMDB 0xC003201FL
#define evtUnhandledExceptionShutdownMDB 0xC0032020L
#define evtUnhandledExceptionMsgSaveChanges 0xC0032021L
#define evtUnhandledExceptionDelete 0xC0032022L
#define evtUnhandledExceptionQuery 0xC0032023L
#define evtUnhandledExceptionFolderSaveChanges 0xC0032024L
#define evtCorruptedCalendar 0xC0032025L
#define evtRebuildCalendar 0x80032026L
#define evtCheckPrimaryCalendar 0x80032027L
#define evtExpandMasterPF 0xC0032028L
#define evtCorruptedPFCalendar 0xC0032029L
#define evtRebuildPFCalendar 0x8003202AL
#define evtMovingMailboxCallbackFailed 0x8003202BL
#endif
