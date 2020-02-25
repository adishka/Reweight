//____________________________________________________________________________
/*
 Copyright (c) 2003-2017, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

          Jim Dobson <J.Dobson07 \at imperial.ac.uk>
          Imperial College London

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TFile.h>
#include <TMath.h>
#include <TNtupleD.h>
#include <cstdlib>
#include <sstream>

// GENIE/Generator includes
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Registry/Registry.h"

// GENIE/Reweight includes
#include "RwCalculators/GReWeightUtils.h"
#include "RwCalculators/GReWeightXSecEmpiricalMEC.h"
#include "RwFramework/GSystSet.h"
#include "RwFramework/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;
using std::ostringstream;

GReWeightXSecEmpiricalMEC::GReWeightXSecEmpiricalMEC() :
GReWeightModel("EmpiricalMEC"),
fManualModelName(),
fManualModelType()
{ this->Init(); }

GReWeightXSecEmpiricalMEC::GReWeightXSecEmpiricalMEC(std::string model, std::string type) :
GReWeightModel("EmpiricalMEC"),
fManualModelName(model),
fManualModelType(type)
{ this->Init(); }

GReWeightXSecEmpiricalMEC::~GReWeightXSecEmpiricalMEC() {
#ifdef _G_REWEIGHT_EmpMEC_DEBUG_
  fTestFile->cd();
  fTestNtp->Write();
  fTestFile->Close();
  delete fTestFile;
#endif
}
//_______________________________________________________________________________________
void GReWeightXSecEmpiricalMEC::Init(void) {
  AlgId id("genie::EmpiricalMECPXSec2015", "Default");

  AlgId twk_id(id);
  if (fManualModelName.size()) {
    twk_id = AlgId(fManualModelName,fManualModelType);
  }

  AlgFactory *algf = AlgFactory::Instance();

  Algorithm *algdef = algf->AdoptAlgorithm(id);
  fXSecModelDef = dynamic_cast<XSecAlgorithmI *>(algdef);
  fXSecModelDef->AdoptSubstructure();

  Algorithm *alg = algf->AdoptAlgorithm(twk_id);
  fXSecModel = dynamic_cast<XSecAlgorithmI *>(alg);
  fXSecModel->AdoptSubstructure();

  fXSecModelConfig = new Registry(fXSecModel->GetConfig());
  fMq2d_Def = fXSecModelDef->GetConfig().GetDouble("EmpiricalMEC-Mq2d");
  fMass_Def = fXSecModelDef->GetConfig().GetDouble("EmpiricalMEC-Mass");
  fWidth_Def = fXSecModelDef->GetConfig().GetDouble("EmpiricalMEC-Width");
  fAPower_Def = fXSecModelDef->GetConfig().GetDouble("EmpiricalMEC-APower");
  fFracCCQE_Def = fXSecModelDef->GetConfig().GetDouble("EmpiricalMEC-FracCCQE");
  fFracPN_NC_Def =
      fXSecModelDef->GetConfig().GetDouble("EmpiricalMEC-FracPN_NC");
  fFracPN_CC_Def =
      fXSecModelDef->GetConfig().GetDouble("EmpiricalMEC-FracPN_CC");
  fFracNCQE_Def = fXSecModelDef->GetConfig().GetDouble("EmpiricalMEC-FracNCQE");
  fFracPN_EM_Def =
      fXSecModelDef->GetConfig().GetDouble("EmpiricalMEC-FracPN_EM");
  fFracEMQE_Def = fXSecModelDef->GetConfig().GetDouble("EmpiricalMEC-FracEMQE");

  fMq2d_TwkDial = 0;
  fMass_TwkDial = 0;
  fWidth_TwkDial = 0;
  fAPower_TwkDial = 0;
  fFracPN_NC_TwkDial = 0;
  fFracPN_CC_TwkDial = 0;
  fFracCCQE_TwkDial = 0;
  fFracNCQE_TwkDial = 0;
  fFracPN_EM_TwkDial = 0;
  fFracEMQE_TwkDial = 0;

#ifdef _G_REWEIGHT_EmpMEC_DEBUG_
  fTestFile = new TFile("./empmec_reweight_test.root", "recreate");
  fTestNtp = new TNtupleD("testntp", "", "E:Q2:wght");
#endif
}
//_______________________________________________________________________________________
bool GReWeightXSecEmpiricalMEC::AppliesTo(ScatteringType_t type,
                                          bool /* is_cc */ ) const {
  /*if (type==kScMEC) { // && is_cc ?
    return true;
  }
  return false;
  */ //adi
  return true;
}
//_______________________________________________________________________________________
bool GReWeightXSecEmpiricalMEC::IsHandled(GSyst_t syst) const {
  if ((syst == kXSecTwkDial_EmpMEC_Mq2d) ||
      (syst == kXSecTwkDial_EmpMEC_Mass) ||
      (syst == kXSecTwkDial_EmpMEC_Width) ||
      (syst == kXSecTwkDial_EmpMEC_APower) ||
      (syst == kXSecTwkDial_EmpMEC_FracPN_NC) ||
      (syst == kXSecTwkDial_EmpMEC_FracPN_CC) ||
      (syst == kXSecTwkDial_EmpMEC_FracCCQE) ||
      (syst == kXSecTwkDial_EmpMEC_FracNCQE) ||
      (syst == kXSecTwkDial_EmpMEC_FracPN_EM) ||
      (syst == kXSecTwkDial_EmpMEC_FracEMQE)) {
    return true;
  }

  return false;
}
//_______________________________________________________________________________________
void GReWeightXSecEmpiricalMEC::SetSystematic(GSyst_t syst, double twk_dial) {
  if (!this->IsHandled(syst)) {
    return;
  }

  switch (syst) {
  case (kXSecTwkDial_EmpMEC_Mq2d): {
    fMq2d_TwkDial = twk_dial;
    break;
  }
  case (kXSecTwkDial_EmpMEC_Mass): {
    fMass_TwkDial = twk_dial;
    break;
  }
  case (kXSecTwkDial_EmpMEC_Width): {
    fWidth_TwkDial = twk_dial;
    break;
  }
  case (kXSecTwkDial_EmpMEC_APower): {
    fAPower_TwkDial = twk_dial;
    break;
  }
  case (kXSecTwkDial_EmpMEC_FracPN_NC): {
    fFracPN_NC_TwkDial = twk_dial;
    break;
  }
  case (kXSecTwkDial_EmpMEC_FracPN_CC): {
    fFracPN_CC_TwkDial = twk_dial;
    break;
  }
  case (kXSecTwkDial_EmpMEC_FracCCQE): {
    fFracCCQE_TwkDial = twk_dial;
    break;
  }
  case (kXSecTwkDial_EmpMEC_FracNCQE): {
    fFracNCQE_TwkDial = twk_dial;
    break;
  }
  case (kXSecTwkDial_EmpMEC_FracPN_EM): {
    fFracPN_EM_TwkDial = twk_dial;
    break;
  }
  case (kXSecTwkDial_EmpMEC_FracEMQE): {
    fFracEMQE_TwkDial = twk_dial;
    break;
  }
  default:
    throw;
  }
}
//_______________________________________________________________________________________
void GReWeightXSecEmpiricalMEC::Reset(void) {
  fMq2d_Curr = fMq2d_Def;
  fMass_Curr = fMass_Def;
  fWidth_Curr = fWidth_Def;
  fAPower_Curr = fAPower_Def;
  fFracPN_NC_Curr = fFracPN_NC_Def;
  fFracPN_CC_Curr = fFracPN_CC_Def;
  fFracCCQE_Curr = fFracCCQE_Def;
  fFracNCQE_Curr = fFracNCQE_Def;
  fFracPN_EM_Curr = fFracPN_EM_Def;
  fFracEMQE_Curr = fFracEMQE_Def;

  this->Reconfigure();
}
//_______________________________________________________________________________________
void GReWeightXSecEmpiricalMEC::Reconfigure(void) {

  bool fMq2d_HasTwk = (fabs(fMq2d_TwkDial) > genie::controls::kASmallNum);
  bool fMass_HasTwk = (fabs(fMass_TwkDial) > genie::controls::kASmallNum);
  bool fWidth_HasTwk = (fabs(fWidth_TwkDial) > genie::controls::kASmallNum);
  bool fAPower_HasTwk = (fabs(fAPower_TwkDial) > genie::controls::kASmallNum);
  bool fFracPN_NC_HasTwk =
      (fabs(fFracPN_NC_TwkDial) > genie::controls::kASmallNum);
  bool fFracPN_CC_HasTwk =
      (fabs(fFracPN_CC_TwkDial) > genie::controls::kASmallNum);
  bool fFracCCQE_HasTwk =
      (fabs(fFracCCQE_TwkDial) > genie::controls::kASmallNum);
  bool fFracNCQE_HasTwk =
      (fabs(fFracNCQE_TwkDial) > genie::controls::kASmallNum);
  bool fFracPN_EM_HasTwk =
      (fabs(fFracPN_EM_TwkDial) > genie::controls::kASmallNum);
  bool fFracEMQE_HasTwk =
      (fabs(fFracEMQE_TwkDial) > genie::controls::kASmallNum);
  fAnyTwk = fMq2d_HasTwk || fMass_HasTwk || fWidth_HasTwk ||
            fFracPN_NC_HasTwk || fFracPN_CC_HasTwk || fFracCCQE_HasTwk || fAPower_HasTwk ||
            fFracNCQE_HasTwk || fFracPN_EM_HasTwk || fFracEMQE_HasTwk;

  if (!fAnyTwk) {
    return;
  }
  Registry r("GReWeightXSecEmpiricalMEC",false);
  //Registry &r = *fXSecModelConfig;
  GSystUncertainty *uncert_inst = GSystUncertainty::Instance();

  if (fMq2d_HasTwk) {
    int sign = utils::rew::Sign(fMq2d_TwkDial);
    double fracerr = uncert_inst->OneSigmaErr(kXSecTwkDial_EmpMEC_Mq2d, sign);
    fMq2d_Curr = fMq2d_Def * (1. + fMq2d_TwkDial * fracerr);
    fMq2d_Curr = TMath::Max(0., fMq2d_Curr); 
    r.Set("EmpiricalMEC-Mq2d", fMq2d_Curr);
  }
  if (fMass_HasTwk) {
    int sign = utils::rew::Sign(fMass_TwkDial);
    double fracerr = uncert_inst->OneSigmaErr(kXSecTwkDial_EmpMEC_Mass, sign);
    fMass_Curr = fMass_Def * (1. + fMass_TwkDial * fracerr);
    fMass_Curr = TMath::Max(0., fMass_Curr); 
    r.Set("EmpiricalMEC-Mass", fMass_Curr);
  }

  if (fWidth_HasTwk) {
    int sign = utils::rew::Sign(fWidth_TwkDial);
    double fracerr = uncert_inst->OneSigmaErr(kXSecTwkDial_EmpMEC_Width, sign);
    fWidth_Curr = fWidth_Def * (1. + fWidth_TwkDial * fracerr);
    fWidth_Curr = TMath::Max(0., fWidth_Curr); 
    r.Set("EmpiricalMEC-Width", fWidth_Curr);
  }
  if (fAPower_HasTwk) {
    int sign = utils::rew::Sign(fAPower_TwkDial);
    double fracerr =
        uncert_inst->OneSigmaErr(kXSecTwkDial_EmpMEC_APower, sign);
    fAPower_Curr = fAPower_Def * (1. + fAPower_TwkDial * fracerr);
    fAPower_Curr = TMath::Max(0., fAPower_Curr);
    fAPower_Curr = TMath::Min(1., fAPower_Curr);
    r.Set("EmpiricalMEC-APower", fAPower_Curr);
  }
  if (fFracPN_NC_HasTwk) {
    int sign = utils::rew::Sign(fFracPN_NC_TwkDial);
    double fracerr =
        uncert_inst->OneSigmaErr(kXSecTwkDial_EmpMEC_FracPN_NC, sign);
    fFracPN_NC_Curr = fFracPN_NC_Def * (1. + fFracPN_NC_TwkDial * fracerr);
    fFracPN_NC_Curr = TMath::Max(0., fFracPN_NC_Curr);
    fFracPN_NC_Curr = TMath::Min(1., fFracPN_NC_Curr);
    r.Set("EmpiricalMEC-FracPN_NC", fFracPN_NC_Curr);
  }
  if (fFracPN_CC_HasTwk) {
    int sign = utils::rew::Sign(fFracPN_CC_TwkDial);
    double fracerr =
        uncert_inst->OneSigmaErr(kXSecTwkDial_EmpMEC_FracPN_CC, sign);
    fFracPN_CC_Curr = fFracPN_CC_Def * (1. + fFracPN_CC_TwkDial * fracerr);
    fFracPN_CC_Curr = TMath::Max(0., fFracPN_CC_Curr);
    fFracPN_CC_Curr = TMath::Min(1., fFracPN_CC_Curr);
    r.Set("EmpiricalMEC-FracPN_CC", fFracPN_CC_Curr);
  }
  if (fFracCCQE_HasTwk) {
    int sign = utils::rew::Sign(fFracCCQE_TwkDial);
    double fracerr =
        uncert_inst->OneSigmaErr(kXSecTwkDial_EmpMEC_FracCCQE, sign);
    fFracCCQE_Curr = fFracCCQE_Def * (1. + fFracCCQE_TwkDial * fracerr);
    fFracCCQE_Curr = TMath::Max(0., fFracCCQE_Curr);
    fFracCCQE_Curr = TMath::Min(1., fFracCCQE_Curr);
    r.Set("EmpiricalMEC-FracCCQE", fFracCCQE_Curr);
  }
  if (fFracNCQE_HasTwk) {
    int sign = utils::rew::Sign(fFracNCQE_TwkDial);
    double fracerr =
        uncert_inst->OneSigmaErr(kXSecTwkDial_EmpMEC_FracNCQE, sign);
    fFracNCQE_Curr = fFracNCQE_Def * (1. + fFracNCQE_TwkDial * fracerr);
    fFracNCQE_Curr = TMath::Max(0., fFracNCQE_Curr);
    fFracNCQE_Curr = TMath::Min(1., fFracNCQE_Curr);
    r.Set("EmpiricalMEC-FracNCQE", fFracNCQE_Curr);
  }
  if (fFracPN_EM_HasTwk) {
    int sign = utils::rew::Sign(fFracPN_EM_TwkDial);
    double fracerr =
        uncert_inst->OneSigmaErr(kXSecTwkDial_EmpMEC_FracPN_EM, sign);
    fFracPN_EM_Curr = fFracPN_EM_Def * (1. + fFracPN_EM_TwkDial * fracerr);
    fFracPN_EM_Curr = TMath::Max(0., fFracPN_EM_Curr);
    fFracPN_EM_Curr = TMath::Min(1., fFracPN_EM_Curr);
    r.Set("EmpiricalMEC-FracPN_EM", fFracPN_EM_Curr);
  }
  if (fFracEMQE_HasTwk) {
    int sign = utils::rew::Sign(fFracEMQE_TwkDial);
    double fracerr =
        uncert_inst->OneSigmaErr(kXSecTwkDial_EmpMEC_FracEMQE, sign);
    fFracEMQE_Curr = fFracEMQE_Def * (1. + fFracEMQE_TwkDial * fracerr);
    fFracEMQE_Curr = TMath::Max(0., fFracEMQE_Curr);
    fFracEMQE_Curr = TMath::Min(1., fFracEMQE_Curr);
    r.Set("EmpiricalMEC-FracEMQE", fFracEMQE_Curr);
  }
  fXSecModel->Configure(r);
}
//_______________________________________________________________________________________
double GReWeightXSecEmpiricalMEC::CalcWeight(const genie::EventRecord &event) {
  
  bool is_mec = event.Summary()->ProcInfo().IsMEC();
  if (!is_mec)
    return 1.0;

  if (!fAnyTwk)
    return 1.0;

  double new_weight = 1.; 
  double weightNorm = 1.;
  double weightPN   = 1.;

  bool tweakedNorm = (TMath::Abs(fFracEMQE_Curr) > controls::kASmallNum);
  if(tweakedNorm){
    weightNorm = CalcWeightNorm(event);
  }
  
  bool tweakedPN = (TMath::Abs(fFracPN_EM_Curr) > controls::kASmallNum);
  if(tweakedPN){
    weightPN = CalcWeightPN(event);
  }

  Interaction *interaction = event.Summary();
  interaction->KinePtr()->UseSelectedKinematics();

  double old_xsec = fXSecModelDef->XSec(interaction, kPSWQ2fE);
  double old_weight = event.Weight();
  double new_xsec = fXSecModel->XSec(interaction, kPSWQ2fE);
  new_weight = old_weight * (new_xsec / old_xsec);
  new_weight *= weightNorm * weightPN;
 
  interaction->KinePtr()->ClearRunningValues();
  return new_weight;
}
//_______________________________________________________________________________________
double GReWeightXSecEmpiricalMEC::CalcWeightNorm(const genie::EventRecord &event) {
  bool tweaked = (TMath::Abs(fFracEMQE_Curr) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  double wght = fFracEMQE_Curr/fFracEMQE_Def;
  return wght;
} 
//_______________________________________________________________________________________
double GReWeightXSecEmpiricalMEC::CalcWeightPN(const genie::EventRecord &event) {
  bool tweaked = (TMath::Abs(fFracPN_EM_Curr) > controls::kASmallNum);
  if(!tweaked) return 1.0;

  Interaction *interaction = event.Summary();
  int nucleon_cluster_pdg = interaction->InitState().Tgt().HitNucPdg();
  double fFracEMQE_cluster_Curr=0.;
  if(nucleon_cluster_pdg==2000000200) fFracEMQE_cluster_Curr = 0.5*(1-fFracPN_EM_Curr);  //n+n
  if(nucleon_cluster_pdg==2000000201) fFracEMQE_cluster_Curr = fFracPN_EM_Curr;  //n+p
  if(nucleon_cluster_pdg==2000000202) fFracEMQE_cluster_Curr = 0.5*(1-fFracPN_EM_Curr);  //p+p
  double fFracEMQE_cluster_Def=0.;
  if(nucleon_cluster_pdg==2000000200) fFracEMQE_cluster_Def = 0.5*(1-fFracPN_EM_Def);  //n+n
  if(nucleon_cluster_pdg==2000000201) fFracEMQE_cluster_Def = fFracPN_EM_Def;  //n+p
  if(nucleon_cluster_pdg==2000000202) fFracEMQE_cluster_Def = 0.5*(1-fFracPN_EM_Def);  //p+p
  double wght = fFracEMQE_cluster_Curr/fFracEMQE_cluster_Def;
  interaction->KinePtr()->ClearRunningValues();
  return wght;
} 
//_______________________________________________________________________________________

