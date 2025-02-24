#include "R3BSofSciRawPosPar.h"

#include "FairLogger.h"
#include "FairParamList.h"

#include "TArrayF.h"
#include "TMath.h"
#include "TString.h"

#include <iostream>

#define MAX_RAWPOSPAR 10

using std::cout;
using std::endl;

// ---- Standard Constructor ---------------------------------------------------
R3BSofSciRawPosPar::R3BSofSciRawPosPar(const char* name, const char* title, const char* context)
    : FairParGenericSet(name, title, context)
    , fNumDets(4)
    , fNumPmts(3)
    , fNumSignals(4)
    , fNumParsPerSignal(2)
{
    fAllSignalsRawPosParams = new TArrayF(fNumSignals * fNumParsPerSignal);
}

// ----  Destructor ------------------------------------------------------------
R3BSofSciRawPosPar::~R3BSofSciRawPosPar()
{
    clear();
    if (fAllSignalsRawPosParams)
    {
        delete fAllSignalsRawPosParams;
    }
}

// ----  Method clear ----------------------------------------------------------
void R3BSofSciRawPosPar::clear()
{
    status = kFALSE;
    resetInputVersions();
}

// ----  Method putParams ------------------------------------------------------
void R3BSofSciRawPosPar::putParams(FairParamList* list)
{
    LOG(INFO) << "R3BSofSciRawPosPar::putParams() called";
    if (!list)
    {
        return;
    }

    Int_t array_size = fNumSignals * fNumParsPerSignal;
    LOG(INFO) << "Array Size: " << array_size;

    fAllSignalsRawPosParams->Set(array_size);

    list->add("RawPosPar", *fAllSignalsRawPosParams);
    list->add("nDetectorsRawPosPar", fNumDets);
    list->add("nChannelsRawPosPar", fNumPmts);
    list->add("nSignalsRawPosPar", fNumSignals);
    list->add("nRawPosParsPerSignal", fNumParsPerSignal);
}

// ----  Method getParams ------------------------------------------------------
Bool_t R3BSofSciRawPosPar::getParams(FairParamList* list)
{
    LOG(INFO) << "R3BSofSciRawPosPar::getParams() called";
    if (!list)
    {
        return kFALSE;
    }
    if (!list->fill("nDetectorsRawPosPar", &fNumDets))
    {
        return kFALSE;
    }
    if (!list->fill("nChannelsRawPosPar", &fNumPmts))
    {
        return kFALSE;
    }
    if (!list->fill("nSignalsRawPosPar", &fNumSignals))
    {
        return kFALSE;
    }
    if (!list->fill("nRawPosParsPerSignal", &fNumParsPerSignal))
    {
        return kFALSE;
    }

    Int_t array_size = fNumSignals * fNumParsPerSignal;
    LOG(INFO) << "Array Size: " << array_size;
    fAllSignalsRawPosParams->Set(array_size);

    if (!(list->fill("RawPosPar", fAllSignalsRawPosParams)))
    {
        LOG(INFO) << "---Could not initialize fAllSignalsRawPosParams";
        return kFALSE;
    }

    return kTRUE;
}

// ----  Method printParams ----------------------------------------------------
void R3BSofSciRawPosPar::printParams()
{
    LOG(INFO) << "R3BSofSciRawPosPar: SofSciRawPos Parameters: ";
    Int_t array_size = fNumSignals * fNumParsPerSignal;

    cout << "--- --------------------------------------------" << endl;
    cout << "--- Single Tcal Parameters :  " << endl;
    cout << "--- --------------------------------------------" << endl;
    for (Int_t param = 0; param < array_size; param++)
    {
        cout << "LIMIT " << param << " = " << fAllSignalsRawPosParams->GetAt(param) << endl;
    }
}
