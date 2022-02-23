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

#include "R3BCalifaP2pAnalysis.h"


// R3BCalifaP2pAnalysis::Constructor
R3BCalifaP2pAnalysis::R3BCalifaP2pAnalysis()
    : FairTask("R3BCalifaP2pAnalysis")
    , fProtonEnergy(20000)
    , fBeta(0.)
    , fPhiSpread(20.)
    , fOnline(kFALSE)
    , fCalifaHitData(NULL)
    , fCalifaP2pData(NULL)
{
}

R3BCalifaP2pAnalysis::~R3BCalifaP2pAnalysis()
{
    R3BLOG(DEBUG1, "");
    if (fCalifaP2pData)
        delete fCalifaP2pData;
}



InitStatus R3BCalifaP2pAnalysis::Init()
{
    R3BLOG(INFO, "Init P2p task");

    // INPUT DATA
    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        R3BLOG(FATAL, "FairRootManager not found");
        return kFATAL;
    }

    fCalifaHitData = (TClonesArray*)rootManager->GetObject("CalifaHitData");
    if (!fCalifaHitData)
    {
        R3BLOG(FATAL, "CalifaHitData not found");
        return kFATAL;
    }

    // OUTPUT DATA
    // Calibrated data
    fCalifaP2pData = new TClonesArray("R3BCalifaP2pData",5);

    rootManager->Register("CalifaP2pData", "CALIFA p2p p2pn Candidates", fCalifaP2pData, !fOnline);

    return kSUCCESS;
}


void R3BCalifaP2pAnalysis::Exec(Option_t* option)
{
    // Reset entries in output arrays, local arrays
    Reset();

    //std::cout<<"TEST"<<std::endl;
    //std::cout<<"P1: "<<energyP1<<std::endl;

    // Reading the Input -- Hit Data --
    nHits = fCalifaHitData->GetEntries();
    if (nHits<=2)
        return;

    R3BCalifaHitData** hitData = new R3BCalifaHitData*[nHits];


    for (Int_t i = 0; i < nHits; i++)
    {
        hitData[i] = (R3BCalifaHitData*)(fCalifaHitData->At(i));
        //Here was some stuff to look for corrupted data. Maybe I do not need this
        tmpProtonPhi = hitData[i]->GetPhi();
        tmpProtonTheta = hitData[i]->GetTheta()*TMath::RadToDeg();
        tmpProtonEnergy = hitData[i]->GetEnergy(); 
        //tmpProtonEnergy=tmpProtonEnergy*dopplerCorrection(tmpProtonTheta,fBeta);
        tmpTime = hitData[i]->GetTime();
       
       // Check if the energy of the hit is above the energy threshold for a real proton hit
       if(tmpProtonEnergy > fProtonEnergy)
       {
           //This if cascade checks, if we have exactly two proton hits. Everything else can not be used!
           if(energyP1>DEFAULTDOUBLE && energyP2>DEFAULTDOUBLE)
           {
               multiEventFlag=1; //There are more then two proton hits -> Can not be used
               break;
           }
           if(energyP1>DEFAULTDOUBLE && energyP2==DEFAULTDOUBLE)
           {
               //There is already one proton hit and this is the second proton hit
               energyP2= tmpProtonEnergy;
               thetaP2=tmpProtonTheta;
               phiP2=tmpProtonPhi;
           }
           if(energyP1==DEFAULTDOUBLE && energyP2==DEFAULTDOUBLE)
           {
               //We have a first proton hit
               energyP1=tmpProtonEnergy;
               thetaP1=tmpProtonTheta;
               phiP1=tmpProtonPhi;
           }  
           if(multiEventFlag==0 && energyP2 != DEFAULTDOUBLE)
                    {
                        //califaP2PCandidate++;
                       //std::cout<<"P2: "<<energyP2<<std::endl;
                       //std::cout<<"T1: "<<thetaP1<<std::endl;
                       //std::cout<<"T2: "<<thetaP2<<std::endl;
                       //std::cout<<"Phi1: "<<phiP1<<std::endl;
                       //std::cout<<"Phi2: "<<phiP2<<std::endl;
                       //std::cout<<"SumEn: "<<sumEnergy<<std::endl;
                        //Tranform the angles from rad to degree
                        //thetaP1=thetaP1;
                        //thetaP2=thetaP2;
                        phiP1=phiP1*TMath::RadToDeg();
                        phiP2=phiP2*TMath::RadToDeg();
                        //Calculate the angles between the two protons and the sum energy
                        thetaP2P=thetaP1+thetaP2;
                        phiP2P=TMath::Abs(phiP1-phiP2);
                        sumEnergy=energyP1+energyP2;
                      
                        if(phiP2P>= 180.-fPhiSpread && phiP2P<= 180.+fPhiSpread)
                        {
                            //std::cout<<"P2: "<<energyP2<<std::endl;
                            //std::cout<<"T1: "<<thetaP1<<std::endl;
                            //std::cout<<"T2: "<<thetaP2<<std::endl;
                            //std::cout<<"Phi1: "<<phiP1<<std::endl;
                            //std::cout<<"Phi2: "hitData[i]->GetTime()<<phiP2<<std::endl;
                            //std::cout<<"SumEn: "<<sumEnergy<<std::endl;
                            //std::cout<<"SumTheta: "<<thetaP2P<<std::endl;
                            //std::cout<<"time: "<<hitData[i]->GetTime()<<std::endl;
                            AddP2pData(energyP1, energyP2, phiP1, phiP2, thetaP1, thetaP2, tmpTime, thetaP2P, phiP2P, sumEnergy);
                            //std::cout<<"ThetaSum: "<<fCalifaP2pData->GetThetaSum()<<std::endl;
                        }

                    }
       }
       //Clear the variables for the next iteration
       
       tmpProtonPhi=DEFAULTDOUBLE;
       tmpProtonTheta=DEFAULTDOUBLE;
       tmpProtonEnergy=DEFAULTDOUBLE; 
       tmpTime=DEFAULTINTEGER; 
        
    }

    if (hitData)
        delete[] hitData;
    return;
}

void R3BCalifaP2pAnalysis::Reset()
{
    nHits=DEFAULTINTEGER;

    energyP2=DEFAULTDOUBLE;
    thetaP2=DEFAULTDOUBLE;
    phiP2=DEFAULTDOUBLE;

    energyP1=DEFAULTDOUBLE;
    thetaP1=DEFAULTDOUBLE;
    phiP1=DEFAULTDOUBLE;

    thetaP2P=DEFAULTDOUBLE;
    phiP2P=DEFAULTDOUBLE;
    sumEnergy=DEFAULTDOUBLE;

    multiEventFlag=0;
    R3BLOG(DEBUG, "Clearing R3BCalifaP2pData Structure");
    if (fCalifaP2pData)
        fCalifaP2pData->Clear();
  
}

R3BCalifaP2pData* R3BCalifaP2pAnalysis::AddP2pData(Double_t energy1,
                                Double_t energy2,
                                Double_t phi1,
                                Double_t phi2,
                                Double_t theta1,
                                Double_t theta2,
                                ULong64_t time,
                                Double_t thetaSum,
                                Double_t phiSum,
                                Double_t energySum)
{
    // It fills the R3BCalifaCrystalCalData
    TClonesArray& clref = *fCalifaP2pData;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BCalifaP2pData(energy1, energy2, energySum, phi1, phi2, phiSum, theta1, theta2, thetaSum, time);
}

Double_t dopplerCorrection(Double_t theta, Double_t BETA)
{
    return (1.-BETA*TMath::Cos(theta))/TMath::Sqrt(1.-BETA*BETA);
}

ClassImp(R3BCalifaP2pAnalysis);
