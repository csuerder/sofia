/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// -----------------------------------------------------------------------------
// -----                                                                   -----
// -----                     R3BCalifaP2pAnalysis                    -----
// -----                Created 16/02/2022 by C. Suerder             -----
// -----                                                                   -----
// -----------------------------------------------------------------------------

#ifndef R3BCALIFAP2PANALYSIS_H
#define R3BCALIFAP2PANALYSIS_H

#include "FairTask.h"

#include "R3BCalifaP2pData.h"

#include "TClonesArray.h"
#include "TMath.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRuntimeDb.h"

#include "R3BCalifaHitData.h"
#include "R3BCalifaP2pData.h"



#include "R3BLogger.h"

class TClonesArray;

class R3BCalifaP2pAnalysis : public FairTask
{
  public:
    /** Default constructor **/
    R3BCalifaP2pAnalysis();

    /** Destructor **/
    virtual ~R3BCalifaP2pAnalysis();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* option);

    /** Virtual method Reset **/
    virtual void Reset();

    //virtual void SetParContainers();

    // Fair specific
    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    //virtual InitStatus ReInit();

    /** Accessor to select online mode **/
    void SetOnline(Bool_t option) { fOnline = option; }
    //void SetThetaSpread(Double_t theta){ fThetaSpread = theta; }
    void SetPhiSpread(Double_t phi){ fPhiSpread = phi; }
    void SetBeta(Double_t beta){ fBeta = beta; }
    void SetMinEnergy( Double_t pEne ){ fProtonEnergy = pEne; }

    //void GetTheta(){ return fTheta; }
    Double_t GetPhiSpread(){ return fPhiSpread; }
    Double_t GetBeta(){ return fBeta; }
    Double_t GetMinEnergy(){ return fProtonEnergy; }


  private:
    //void SetParameter();
    Int_t DEFAULTINTEGER=0;
    Double_t DEFAULTDOUBLE=0.;

    // Don't store data for online
    Bool_t fOnline;
    //Double_t fThetaSpread;
    
    Double_t fPhiSpread;
    Double_t fBeta;                 //----------------> Get it from FrsData?
    Double_t fProtonEnergy;

    TClonesArray* fCalifaHitData; /**< Array with CALIFA Mapped- input data. >*/
    TClonesArray* fCalifaP2pData; /**< Array with CALIFA Cal- output data. >*/

    Int_t nHits=DEFAULTINTEGER;
    Double_t tmpProtonPhi=DEFAULTDOUBLE;        //Temporary store the phi of the proton hit
    Double_t tmpProtonTheta=DEFAULTDOUBLE;      //Temporary store the theta of the proton hit
    Double_t tmpProtonEnergy=DEFAULTDOUBLE;     //Temporary store the doppler corrected energy of the proton hit
    Int_t multiEventFlag=DEFAULTDOUBLE;                     //If 1, we have more then two proton hits and the event is discarded
    Long64_t tmpTime=DEFAULTINTEGER;

    Double_t energyP1=DEFAULTDOUBLE;    //Energy of proton 1 in the rest frame
    Double_t thetaP1=DEFAULTDOUBLE;     //Theta of proton 1
    Double_t phiP1=DEFAULTDOUBLE;       //Phi of proton 1

    Double_t energyP2=DEFAULTDOUBLE;    //Energy of proton 2 in the rest frame
    Double_t thetaP2=DEFAULTDOUBLE;     //Theta of proton 2
    Double_t phiP2=DEFAULTDOUBLE;       //Phi of proton 2

    Double_t thetaP2P=DEFAULTDOUBLE;            //The theta between proton 1 and 2 (only true, if phi is ~180 degree)
    Double_t phiP2P=DEFAULTDOUBLE;              //The phi between proton 1 and 2
    Double_t sumEnergy=DEFAULTDOUBLE; 

    R3BCalifaP2pData* AddP2pData(Double_t fEnergy1,
                                Double_t fEnergy2,
                                Double_t fPhi1,
                                Double_t fPhi2,
                                Double_t fTheta1,
                                Double_t fTheta2,
                                ULong64_t fTime,
                                Double_t fThetaSum,
                                Double_t fPhiSum,
                                Double_t fEnergySum);

    Double_t dopplerCorrection(Double_t theta, Double_t BETA);


  public:
    // Class definition
    ClassDef(R3BCalifaP2pAnalysis, 1)
};

#endif /* R3BCalifaP2pAnalysis_H */
