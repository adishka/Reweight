//____________________________________________________________________________
/*!

\class    genie::IBDPrimaryLeptonGenerator

\brief    Generates the final state primary lepton in v IBD interactions.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Corey Reed <cjreed \at nikhef.nl> - October 29, 2009
          using code from the QELKinematicGenerator written by
          Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

\created  October 29, 2009

\cpright  Copyright (c) 2003-2013, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _IBD_PRIMARY_LEPTON_GENERATOR_H_
#define _IBD_PRIMARY_LEPTON_GENERATOR_H_

#include "EVGModules/PrimaryLeptonGenerator.h"

namespace genie {

class IBDPrimaryLeptonGenerator : public PrimaryLeptonGenerator {

public :
  IBDPrimaryLeptonGenerator();
  IBDPrimaryLeptonGenerator(string config);
  virtual ~IBDPrimaryLeptonGenerator();

  // implement the EventRecordVisitorI interface
  void ProcessEventRecord(GHepRecord * event_rec) const;
};

}      // genie namespace
#endif // _IBD_PRIMARY_LEPTON_GENERATOR_H_
