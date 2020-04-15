// ---------------------------------------------------------------------------
// -----                                                                 -----
// -----                      R3BSofFrsData                              -----
// -----                  Created 28/01/2020 by J.L. Rodriguez           -----
// -----                                                                 -----
// ---------------------------------------------------------------------------

#ifndef R3BSofFrsData_H
#define R3BSofFrsData_H
#include "TObject.h"

class R3BSofFrsData : public TObject
{

  public:
    // Default Constructor
    R3BSofFrsData();

    /** Standard Constructor
     *@param fZ      Z of fragments
     *@param fAq     A/q of fragments
     *@param fBeta   Beta of fragments
     *@param fBrho   Brho of fragments
     *@param cs2     Position at S2
     *@param cs8     Position at S8
     *@param xc      Position at Cave
     **/
    R3BSofFrsData(Double_t z, Double_t aq, Double_t beta, Double_t brho, Double_t xs2 = 0., Double_t xc = 0.);

    // Destructor
    virtual ~R3BSofFrsData() {}

    // Getters
    inline const Double_t GetZ() const { return fZ; }
    inline const Double_t GetAq() const { return fAq; }
    inline const Double_t GetBeta() const { return fBeta; }
    inline const Double_t GetBrho() const { return fBrho; }
    inline const Double_t GetXS2() const { return fxs2; }
    //inline const Double_t GetXS8() const { return fxs8; }
    inline const Double_t GetXCave() const { return fxc; }

  protected:
    Double_t fZ, fAq; // ID
    Double_t fBeta, fBrho, fxs2, fxc; //fxs8

  public:
    ClassDef(R3BSofFrsData, 1)
};

#endif
