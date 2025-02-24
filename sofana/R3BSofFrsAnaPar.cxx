// ------------------------------------------------------------------
// -----         R3BSofFrsAnaPar source file                    -----
// -----         Created 27/01/20  by J.L. Rodriguez-Sanchez    -----
// -----         Revised 07/08/20  by R. Taniuchi               -----
// ------------------------------------------------------------------

#include "R3BSofFrsAnaPar.h"
#define MAX_TOFNUM 5

// ---- Standard Constructor ---------------------------------------------------
R3BSofFrsAnaPar::R3BSofFrsAnaPar(const TString& name, const TString& title, const TString& context)
    : FairParGenericSet(name, title, context)
    , fBrho0(0)
    , fNumTof(0)
    , fS2PosCoef(0)
    , fS2PosOffset(0)
    , fNumBrhoCorrPar(0)
{
    fStaSciId = new TArrayI(MAX_TOFNUM);
    fStoSciId = new TArrayI(MAX_TOFNUM);
    fPathLength = new TArrayF(MAX_TOFNUM);
    fTofOffset = new TArrayF(MAX_TOFNUM);
    fUseS2x = new TArrayI(MAX_TOFNUM);
    fBrhoCorrPar = new TArrayF(4);
}

// ----  Destructor ------------------------------------------------------------
R3BSofFrsAnaPar::~R3BSofFrsAnaPar()
{
    clear();
    if (fStaSciId)
        delete fStaSciId;
    if (fStoSciId)
        delete fStoSciId;
    if (fPathLength)
        delete fPathLength;
    if (fTofOffset)
        delete fTofOffset;
    if (fUseS2x)
        delete fUseS2x;
    if (fBrhoCorrPar)
        delete fBrhoCorrPar;
}

// ----  Method clear ----------------------------------------------------------
void R3BSofFrsAnaPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BSofFrsAnaPar::putParams(FairParamList* list)
{
    LOG(INFO) << "R3BSofFrsAnaPar::putParams() called";
    if (!list)
    {
        return;
    }
    list->add("BrhoS2S8", fBrho0);
    list->add("NumTof", fNumTof);

    fStaSciId->Set(fNumTof);
    fStoSciId->Set(fNumTof);
    fPathLength->Set(fNumTof);
    fTofOffset->Set(fNumTof);
    fUseS2x->Set(fNumTof);
    list->add("StaSciId", *fStaSciId);
    list->add("StoSciId", *fStoSciId);
    list->add("PathLength", *fPathLength);
    list->add("TofOffset", *fTofOffset);
    list->add("UseS2x", *fUseS2x);

    list->add("S2PosCoef", fS2PosCoef);
    list->add("S2PosOffset", fS2PosOffset);

    list->add("NumBrhoCorrPar", fNumBrhoCorrPar);
    list->add("BrhoCorrPar", *fBrhoCorrPar);
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BSofFrsAnaPar::getParams(FairParamList* list)
{
    LOG(INFO) << "R3BSofFrsAnaPar::getParams() called";
    if (!list)
    {
        return kFALSE;
    }

    if (!list->fill("BrhoS2S8", &fBrho0))
        return kFALSE;
    if (!list->fill("NumTof", &fNumTof))
        return kFALSE;

    fStaSciId->Set(fNumTof);
    fStoSciId->Set(fNumTof);
    fPathLength->Set(fNumTof);
    fTofOffset->Set(fNumTof);
    fUseS2x->Set(fNumTof);
    if (!list->fill("StaSciId", fStaSciId))
        return kFALSE;
    if (!list->fill("StoSciId", fStoSciId))
        return kFALSE;
    if (!list->fill("PathLength", fPathLength))
        return kFALSE;
    if (!list->fill("TofOffset", fTofOffset))
        return kFALSE;
    if (!list->fill("UseS2x", fUseS2x))
        return kFALSE;

    if (!list->fill("S2PosCoef", &fS2PosCoef))
        return kFALSE;
    if (!list->fill("S2PosOffset", &fS2PosOffset))
        return kFALSE;

    if (!list->fill("NumBrhoCorrPar", &fNumBrhoCorrPar))
        return kFALSE;
    if (!list->fill("BrhoCorrPar", fBrhoCorrPar))
        return kFALSE;

    //printParams();

    return kTRUE;
}

// ----  Method printParams ----------------------------------------------------
void R3BSofFrsAnaPar::printParams()
{
    LOG(INFO) << "R3BSofFrsAnaPar: Frs Analysis Parameters: ";
    LOG(INFO) << "BrhoS2S8: " << fBrho0;
    LOG(INFO) << "NumTof: " << fNumTof;
    for (Int_t i = 0; i < fNumTof; i++)
    {
        LOG(INFO) << "StaSciId(" << i << "): " << fStaSciId->GetAt(i);
        LOG(INFO) << "StoSciId(" << i << "): " << fStoSciId->GetAt(i);
        LOG(INFO) << "PathLength(" << i << "): " << fPathLength->GetAt(i);
        LOG(INFO) << "TofOffset(" << i << "): " << fTofOffset->GetAt(i);
        LOG(INFO) << "UseS2x(" << i << "): " << fUseS2x->GetAt(i);
    }
    LOG(INFO) << "S2PosCoef: " << fS2PosCoef;
    LOG(INFO) << "S2PosOffset: " << fS2PosOffset;
    for (Int_t i = 0; i < fNumBrhoCorrPar; i++)
      {
	LOG(INFO) << "BrhoCorrPar(" << i << "): " << fBrhoCorrPar->GetAt(i);
      }
}
