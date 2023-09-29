// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME GWADict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "src/channel.h"
#include "src/polesearcher.h"
#include "src/bottomonium.h"
#include "src/amplitude.h"
#include "src/filereader.h"
#include "src/observable.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace {
  void TriggerDictionaryInitialization_GWADict_Impl() {
    static const char* headers[] = {
"src/channel.h",
"src/polesearcher.h",
"src/bottomonium.h",
"src/amplitude.h",
"src/filereader.h",
"src/observable.h",
nullptr
    };
    static const char* includePaths[] = {
"/geode2/home/u080/smithwya/Carbonate/root_install/include/",
"/geode2/home/u080/smithwya/Carbonate/GWA-PWA/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "GWADict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "GWADict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "src/channel.h"
#include "src/polesearcher.h"
#include "src/bottomonium.h"
#include "src/amplitude.h"
#include "src/filereader.h"
#include "src/observable.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("GWADict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_GWADict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_GWADict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_GWADict() {
  TriggerDictionaryInitialization_GWADict_Impl();
}
