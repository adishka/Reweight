//____________________________________________________________________________
/*!

\class    genie::rew::GReWeightXSecEMQE

\brief    Reweighting vector form factors in GENIE CCQE neutrino cross
          section calculations.

\author   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

\created  May 24, 2010

\cpright  Copyright (c) 2003-2018, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _G_REWEIGHT_EMQE_H_
#define _G_REWEIGHT_EMQE_H_

//#define _G_REWEIGHT_CCQE_VEC_DEBUG_

#include <map>
#include <string>

// GENIE/Reweight includes
#include "RwCalculators/GReWeightModel.h"


class TFile;
class TNtupleD;

namespace genie {

class XSecAlgorithmI;

namespace rew   {

 class GReWeightXSecEMQE : public GReWeightModel
 {
 public:
   GReWeightXSecEMQE();
  ~GReWeightXSecEMQE();

   // implement the GReWeightI interface
   bool   AppliesTo      (ScatteringType_t type, bool is_cc) const;
   bool   IsHandled      (GSyst_t syst) const;
   void   SetSystematic  (GSyst_t syst, double val);
   void   Reset          (void);
   void   Reconfigure    (void);
   double CalcWeight     (const EventRecord & event);

   // various config options
   void RewE     (bool tf ) { fRewE    = tf;   }
   void RewEbar  (bool tf ) { fRewEbar = tf;   }
   void SetELMvPath   (string p) { fELMvPath     = p;    }
 private:

   void Init (void);

   XSecAlgorithmI * fXSecModel_bba;  ///< EMQE model with BBA05  f/f (default)
   XSecAlgorithmI * fXSecModel_dpl;  ///< EMQE model with dipole f/f ("maximally" tweaked)
   Registry *       fXSecModelConfig; ///< config in tweaked model

   double fFFTwkDial;    ///< tweaking dial (0: bba/default, +1: dipole)

   bool fRewE;
   bool fRewEbar;
   string fELMvPath;       ///<
   double fELMvDef;        ///<
   double fELMvCurr;       ///<

#ifdef _G_REWEIGHT_EMQE_DEBUG_
   TFile *    fTestFile;
   TNtupleD * fTestNtp;
#endif
 };

} // rew   namespace
} // genie namespace

#endif
