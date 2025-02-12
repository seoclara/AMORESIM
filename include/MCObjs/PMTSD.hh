#ifndef PMTSD_HH
#define PMTSD_HH

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// PMTSD                                                               //
//                                                                      //
// Description of the event parameters                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class PMTSD : public TObject {

  private:
    Int_t Nhits;
    Int_t NhitPmts;
    Int_t NtotPmts;

  public:
    PMTSD();
    PMTSD(const PMTSD &orig);
    virtual ~PMTSD() { Clear(); }
    PMTSD &operator=(const PMTSD &orig);

    void Clear(Option_t *option = "");
    Int_t GetNhits() const { return Nhits; }
    Int_t GetNhitPmts() const { return NhitPmts; }
    Int_t GetNtotPmts() const { return NtotPmts; }

    void SetNhits(Int_t nhit) { Nhits = nhit; }
    void SetNhitPmts(Int_t npmt) { NhitPmts = npmt; }
    void SetNtotPmts(Int_t npmt) { NtotPmts = npmt; }

    ClassDef(PMTSD, 2) // Track structure
};

#endif
