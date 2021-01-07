// ------------------------------------------------------------
// -----                  R3BSofFrsFillTree               -----
// -----              Fill SOFIA and FRS tree             ----
// -----           Created 07/08/20 by R. Taniuchi        -----
// ------------------------------------------------------------

#ifndef R3BSofFrsFillTree_H
#define R3BSofFrsFillTree_H

#include "FairTask.h"
#include "TMath.h"
#include "TClonesArray.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

#include "R3BEventHeader.h"
#include "R3BMusicCalData.h"
#include "R3BMusicHitData.h"
#include "R3BSofFrsData.h"
#include "R3BSofMwpcCalData.h"
#include "R3BSofSciCalData.h"
#include "R3BSofSciMappedData.h"
#include "R3BSofSciSingleTcalData.h"
#include "R3BSofSciTcalData.h"
#include "R3BSofTwimCalData.h"
#include "R3BSofTwimHitData.h"
#include "R3BSofTwimHitPar.h"

class TClonesArray;
class R3BEventHeader;

/**
 * This taks reads SCI data and plots online histograms
 */
class R3BSofFrsFillTree : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BSofFrsFillTree();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BSofFrsFillTree(const char* name, Int_t iVerbose = 1);
    
    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BSofFrsFillTree();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t* option);

    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    virtual void FinishEvent();

    /**
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    virtual void FinishTask();

    /**
     * Methods to clean histograms.
     */
    virtual void Reset_Histo();

    /** Virtual method Reset **/
    virtual void Reset();

    void SetNbDetectors(UChar_t ndets) {fNbDetectors = ndets;}
    void SetNbChannels(UChar_t nchs) {fNbChannels = nchs;}
    void SetIdS2(UChar_t id) {fIdS2 = id;}
    void SetIdS8(UChar_t id) {fIdS8 = id;}
    void SetIdCave(UChar_t id) {fIdCave = id;}

    UChar_t GetNbDetectors() {return fNbDetectors;}
    UChar_t GetNbChannels() {return fNbChannels;}
    UChar_t GetIdS2() {return fIdS2;}
    UChar_t GetIdS8() {return fIdS8;}
    UChar_t GetIdCave() {return fIdCave;}

  private:
    TClonesArray* fMappedItemsSci;     /**< Array with mapped items. */
    TClonesArray* fTcalItemsSci;       /**< Array with tcal items. */
    TClonesArray* fSingleTcalItemsSci; /**< Array with tcal items. */
    TClonesArray* fMusHitItems;        /**< Array with MUSIC Hit items. */
    TClonesArray* fMusCalItems;        /**< Array with MUSIC Cal items. */
    TClonesArray* fTwimHitItems;        /**< Array with Twim Hit items. */
    TClonesArray* fCalItemsMwpc0;      /**< Array with cal items of mwpc0. */
    TClonesArray* fFrsData;

    //Scintillators
    UChar_t fNbDetectors;
    UChar_t fNbChannels;
    UChar_t fIdS2;
    UChar_t fIdS8;
    UChar_t fIdCave;

    //Musics
    UChar_t fNumSec;
    UChar_t fNumAnodes;
    UChar_t fNumParams;
    Float_t fZ0 = 0., fZ1 = 0. , fZ2 = 0.; // CalibPar for R3BMUSIC
    TArrayF* CalZParams;
    Float_t fTwimZ0 = 0., fTwimZ1 = 0., fTwimZ2 = 0.; // CalibPar for Twim
    TArrayF* TwimCalZParams;

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header.      */
    UInt_t fNEvents;         /**< Event counter.     */

    // Tree and branches
    TTree* FrsTree;
    Float_t MusicZ = -5000., MusicE = -5000.;
    Float_t MusicDT = -5000.;
    Float_t TwimE = -5000., TwimZ = -5000.;
    UChar_t* multMapSci;
    Float_t* iRawTimeNs;
    Float_t xs2 = -5000., xpos[3] = {-5000.};
    Float_t MusicZ_S2_Cave = -5000., Tof_wTref_S2_Cave = -5000., Beta_S2_Cave = -5000., Brho_S2_Cave = -5000.;
    Float_t MusicZ_S2_S8 = -5000., Tof_wTref_S2_S8 = -5000., Beta_S2_S8 = -5000., Brho_S2_S8 = -5000.;
    Float_t MusicZ_S8_Cave = -5000., Tof_wTref_S8_Cave = -5000., Beta_S8_Cave = -5000., Brho_S8_Cave = -5000.;
    Float_t AoQ_S2_Cave = -5000., AoQ_S2_S8 = -5000., AoQ_S8_Cave = -5000.;
    Float_t TheBeta = -5000., TheBrho = -5000., TheAoQ = -5000.;
    UInt_t tpat = 0, trigger = 0;
    
  public:
    ClassDef(R3BSofFrsFillTree, 1)
};

#endif
