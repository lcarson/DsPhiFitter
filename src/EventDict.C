// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIhomedIhadavizadehdIBc_AnalysisdIDataStrippingdIB_PhiDdIGit_FitdITestdIDsPhiFitterdIsrcdIEventDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
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

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/home/hadavizadeh/Bc_Analysis/DataStripping/B_PhiD/Git_Fit/Test/DsPhiFitter/src/RooHILLdini.C"
#include "/home/hadavizadeh/Bc_Analysis/DataStripping/B_PhiD/Git_Fit/Test/DsPhiFitter/src/RooHORNSdini.C"
#include "/home/hadavizadeh/Bc_Analysis/DataStripping/B_PhiD/Git_Fit/Test/DsPhiFitter/src/RooLITTLEHORNSdini.C"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_RooHILLdini(void *p);
   static void deleteArray_RooHILLdini(void *p);
   static void destruct_RooHILLdini(void *p);
   static void streamer_RooHILLdini(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooHILLdini*)
   {
      ::RooHILLdini *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooHILLdini >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooHILLdini", ::RooHILLdini::Class_Version(), "src/RooHILLdini.h", 27,
                  typeid(::RooHILLdini), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RooHILLdini::Dictionary, isa_proxy, 16,
                  sizeof(::RooHILLdini) );
      instance.SetDelete(&delete_RooHILLdini);
      instance.SetDeleteArray(&deleteArray_RooHILLdini);
      instance.SetDestructor(&destruct_RooHILLdini);
      instance.SetStreamerFunc(&streamer_RooHILLdini);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooHILLdini*)
   {
      return GenerateInitInstanceLocal((::RooHILLdini*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::RooHILLdini*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_RooHORNSdini(void *p);
   static void deleteArray_RooHORNSdini(void *p);
   static void destruct_RooHORNSdini(void *p);
   static void streamer_RooHORNSdini(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooHORNSdini*)
   {
      ::RooHORNSdini *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooHORNSdini >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooHORNSdini", ::RooHORNSdini::Class_Version(), "src/RooHORNSdini.h", 26,
                  typeid(::RooHORNSdini), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RooHORNSdini::Dictionary, isa_proxy, 16,
                  sizeof(::RooHORNSdini) );
      instance.SetDelete(&delete_RooHORNSdini);
      instance.SetDeleteArray(&deleteArray_RooHORNSdini);
      instance.SetDestructor(&destruct_RooHORNSdini);
      instance.SetStreamerFunc(&streamer_RooHORNSdini);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooHORNSdini*)
   {
      return GenerateInitInstanceLocal((::RooHORNSdini*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::RooHORNSdini*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_RooLITTLEHORNSdini(void *p);
   static void deleteArray_RooLITTLEHORNSdini(void *p);
   static void destruct_RooLITTLEHORNSdini(void *p);
   static void streamer_RooLITTLEHORNSdini(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooLITTLEHORNSdini*)
   {
      ::RooLITTLEHORNSdini *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooLITTLEHORNSdini >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooLITTLEHORNSdini", ::RooLITTLEHORNSdini::Class_Version(), "src/RooLITTLEHORNSdini.h", 27,
                  typeid(::RooLITTLEHORNSdini), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RooLITTLEHORNSdini::Dictionary, isa_proxy, 16,
                  sizeof(::RooLITTLEHORNSdini) );
      instance.SetDelete(&delete_RooLITTLEHORNSdini);
      instance.SetDeleteArray(&deleteArray_RooLITTLEHORNSdini);
      instance.SetDestructor(&destruct_RooLITTLEHORNSdini);
      instance.SetStreamerFunc(&streamer_RooLITTLEHORNSdini);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooLITTLEHORNSdini*)
   {
      return GenerateInitInstanceLocal((::RooLITTLEHORNSdini*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::RooLITTLEHORNSdini*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr RooHILLdini::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooHILLdini::Class_Name()
{
   return "RooHILLdini";
}

//______________________________________________________________________________
const char *RooHILLdini::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooHILLdini*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooHILLdini::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooHILLdini*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooHILLdini::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooHILLdini*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooHILLdini::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooHILLdini*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RooHORNSdini::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooHORNSdini::Class_Name()
{
   return "RooHORNSdini";
}

//______________________________________________________________________________
const char *RooHORNSdini::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooHORNSdini*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooHORNSdini::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooHORNSdini*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooHORNSdini::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooHORNSdini*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooHORNSdini::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooHORNSdini*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RooLITTLEHORNSdini::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooLITTLEHORNSdini::Class_Name()
{
   return "RooLITTLEHORNSdini";
}

//______________________________________________________________________________
const char *RooLITTLEHORNSdini::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooLITTLEHORNSdini*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooLITTLEHORNSdini::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooLITTLEHORNSdini*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooLITTLEHORNSdini::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooLITTLEHORNSdini*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooLITTLEHORNSdini::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooLITTLEHORNSdini*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void RooHILLdini::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooHILLdini.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      m.Streamer(R__b);
      a.Streamer(R__b);
      b.Streamer(R__b);
      csi.Streamer(R__b);
      shift.Streamer(R__b);
      sigma.Streamer(R__b);
      ratio_sigma.Streamer(R__b);
      fraction_sigma.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, RooHILLdini::IsA());
   } else {
      R__c = R__b.WriteVersion(RooHILLdini::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      m.Streamer(R__b);
      a.Streamer(R__b);
      b.Streamer(R__b);
      csi.Streamer(R__b);
      shift.Streamer(R__b);
      sigma.Streamer(R__b);
      ratio_sigma.Streamer(R__b);
      fraction_sigma.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_RooHILLdini(void *p) {
      delete ((::RooHILLdini*)p);
   }
   static void deleteArray_RooHILLdini(void *p) {
      delete [] ((::RooHILLdini*)p);
   }
   static void destruct_RooHILLdini(void *p) {
      typedef ::RooHILLdini current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RooHILLdini(TBuffer &buf, void *obj) {
      ((::RooHILLdini*)obj)->::RooHILLdini::Streamer(buf);
   }
} // end of namespace ROOT for class ::RooHILLdini

//______________________________________________________________________________
void RooHORNSdini::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooHORNSdini.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      m.Streamer(R__b);
      a.Streamer(R__b);
      b.Streamer(R__b);
      csi.Streamer(R__b);
      shift.Streamer(R__b);
      sigma.Streamer(R__b);
      ratio_sigma.Streamer(R__b);
      fraction_sigma.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, RooHORNSdini::IsA());
   } else {
      R__c = R__b.WriteVersion(RooHORNSdini::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      m.Streamer(R__b);
      a.Streamer(R__b);
      b.Streamer(R__b);
      csi.Streamer(R__b);
      shift.Streamer(R__b);
      sigma.Streamer(R__b);
      ratio_sigma.Streamer(R__b);
      fraction_sigma.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_RooHORNSdini(void *p) {
      delete ((::RooHORNSdini*)p);
   }
   static void deleteArray_RooHORNSdini(void *p) {
      delete [] ((::RooHORNSdini*)p);
   }
   static void destruct_RooHORNSdini(void *p) {
      typedef ::RooHORNSdini current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RooHORNSdini(TBuffer &buf, void *obj) {
      ((::RooHORNSdini*)obj)->::RooHORNSdini::Streamer(buf);
   }
} // end of namespace ROOT for class ::RooHORNSdini

//______________________________________________________________________________
void RooLITTLEHORNSdini::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooLITTLEHORNSdini.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      m.Streamer(R__b);
      a.Streamer(R__b);
      b.Streamer(R__b);
      csi.Streamer(R__b);
      shift.Streamer(R__b);
      sigma.Streamer(R__b);
      ratio_sigma.Streamer(R__b);
      fraction_sigma.Streamer(R__b);
      shiftg.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, RooLITTLEHORNSdini::IsA());
   } else {
      R__c = R__b.WriteVersion(RooLITTLEHORNSdini::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      m.Streamer(R__b);
      a.Streamer(R__b);
      b.Streamer(R__b);
      csi.Streamer(R__b);
      shift.Streamer(R__b);
      sigma.Streamer(R__b);
      ratio_sigma.Streamer(R__b);
      fraction_sigma.Streamer(R__b);
      shiftg.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_RooLITTLEHORNSdini(void *p) {
      delete ((::RooLITTLEHORNSdini*)p);
   }
   static void deleteArray_RooLITTLEHORNSdini(void *p) {
      delete [] ((::RooLITTLEHORNSdini*)p);
   }
   static void destruct_RooLITTLEHORNSdini(void *p) {
      typedef ::RooLITTLEHORNSdini current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RooLITTLEHORNSdini(TBuffer &buf, void *obj) {
      ((::RooLITTLEHORNSdini*)obj)->::RooLITTLEHORNSdini::Streamer(buf);
   }
} // end of namespace ROOT for class ::RooLITTLEHORNSdini

namespace {
  void TriggerDictionaryInitialization_EventDict_Impl() {
    static const char* headers[] = {
"src/RooHILLdini.C",
"src/RooHORNSdini.C",
"src/RooLITTLEHORNSdini.C",
0
    };
    static const char* includePaths[] = {
"/cvmfs/lhcb.cern.ch/lib/lcg/releases/ROOT/6.06.02-6cc9c/x86_64-slc6-gcc49-opt/include",
"/home/hadavizadeh/Bc_Analysis/DataStripping/B_PhiD/Git_Fit/Test/DsPhiFitter/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "EventDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(RooHILLdini function PDF)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$src/RooHILLdini.C")))  RooHILLdini;
class __attribute__((annotate(R"ATTRDUMP(RooHORNSdini function PDF)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$src/RooHORNSdini.C")))  RooHORNSdini;
class __attribute__((annotate(R"ATTRDUMP(RooLITTLEHORNSdini function PDF)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$src/RooLITTLEHORNSdini.C")))  RooLITTLEHORNSdini;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "EventDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "src/RooHILLdini.C"
#include "src/RooHORNSdini.C"
#include "src/RooLITTLEHORNSdini.C"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"RooHILLdini", payloadCode, "@",
"RooHORNSdini", payloadCode, "@",
"RooLITTLEHORNSdini", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("EventDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_EventDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_EventDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_EventDict() {
  TriggerDictionaryInitialization_EventDict_Impl();
}
