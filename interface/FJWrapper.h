#ifndef FJWrapper_H
#define FJWrapper_H

#include <vector>
#include <TString.h>
#include "FJ_includes.h"
#include "FJJetShape.h"

class FJWrapper
{
 public:
  FJWrapper(const char *name, const char *title);
  virtual ~FJWrapper();

  virtual void  AddInputVector (Double_t px, Double_t py, Double_t pz, Double_t E, Int_t index = -99999);
  virtual void  AddInputVector (const fastjet::PseudoJet& vec,                Int_t index = -99999);
  virtual void  AddInputVectors(const std::vector<fastjet::PseudoJet>& vecs,  Int_t offsetIndex = -99999);
  virtual void  AddInputGhost  (Double_t px, Double_t py, Double_t pz, Double_t E, Int_t index = -99999);
  virtual const char *ClassName()                            const { return "FJWrapper";              }
  virtual void  Clear(const Option_t* /*opt*/ = "");
  virtual void  CopySettingsFrom (const FJWrapper& wrapper);
  virtual void  GetMedianAndSigma(Double_t& median, Double_t& sigma, Int_t remove = 0) const;
  fastjet::ClusterSequenceArea*           GetClusterSequence() const   { return fClustSeq;                 }
  fastjet::ClusterSequence*               GetClusterSequenceSA() const { return fClustSeqSA;               }
  fastjet::ClusterSequenceActiveAreaExplicitGhosts* GetClusterSequenceGhosts() const { return fClustSeqActGhosts; }
  const std::vector<fastjet::PseudoJet>&  GetInputVectors()    const { return fInputVectors;               }
  const std::vector<fastjet::PseudoJet>&  GetInputGhosts()     const { return fInputGhosts;                }
  const std::vector<fastjet::PseudoJet>&  GetInclusiveJets()   const { return fInclusiveJets;              }
  const std::vector<fastjet::PseudoJet>&  GetFilteredJets()    const { return fFilteredJets;               }
  std::vector<fastjet::PseudoJet>         GetJetConstituents(UInt_t idx) const;
  std::vector<fastjet::PseudoJet>         GetFilteredJetConstituents(UInt_t idx) const;
  Double_t                                GetMedianUsedForBgSubtraction() const { return fMedUsedForBgSub; }
  const char*                             GetName()            const { return fName;                       }
  const char*                             GetTitle()           const { return fTitle;                      }
  Double_t                                GetJetArea         (UInt_t idx) const;
  fastjet::PseudoJet                      GetJetAreaVector   (UInt_t idx) const;
  Double_t                                GetFilteredJetArea (UInt_t idx) const;
  fastjet::PseudoJet                      GetFilteredJetAreaVector(UInt_t idx) const;
  Double_t                                GetJetSubtractedPt (UInt_t idx) const;
  virtual std::vector<double>             GetSubtractedJetsPts(Double_t median_pt = -1, Bool_t sorted = kFALSE);
  Bool_t                                  GetLegacyMode()            { return fLegacyMode; }
  Bool_t                                  GetDoFilterArea()          { return fDoFilterArea; }

  const std::vector<fastjet::contrib::GenericSubtractorInfo> GetGenSubtractorInfoJetMass()        const {return fGenSubtractorInfoJetMass        ; }
  const std::vector<fastjet::contrib::GenericSubtractorInfo> GetGenSubtractorInfoJetAngularity()  const {return fGenSubtractorInfoJetAngularity  ; } 
  const std::vector<fastjet::contrib::GenericSubtractorInfo> GetGenSubtractorInfoJetpTD()         const {return fGenSubtractorInfoJetpTD         ; }
  const std::vector<fastjet::contrib::GenericSubtractorInfo> GetGenSubtractorInfoJetCircularity() const {return fGenSubtractorInfoJetCircularity ; }
  const std::vector<fastjet::contrib::GenericSubtractorInfo> GetGenSubtractorInfoJetSigma2()      const {return fGenSubtractorInfoJetSigma2      ; }
  const std::vector<fastjet::contrib::GenericSubtractorInfo> GetGenSubtractorInfoJetConstituent() const {return fGenSubtractorInfoJetConstituent ; } 
  const std::vector<fastjet::contrib::GenericSubtractorInfo> GetGenSubtractorInfoJetLeSub()       const {return fGenSubtractorInfoJetLeSub       ; }
  const std::vector<fastjet::PseudoJet>                      GetConstituentSubtrJets()            const {return fConstituentSubtrJets            ; }
 
  virtual std::vector<double>                                GetGRNumerator()                     const { return fGRNumerator                    ; }
  virtual std::vector<double>                                GetGRDenominator()                   const { return fGRDenominator                  ; }
  virtual std::vector<double>                                GetGRNumeratorSub()                  const { return fGRNumeratorSub                 ; }
  virtual std::vector<double>                                GetGRDenominatorSub()                const { return fGRDenominatorSub               ; }

  virtual Int_t Run();
  virtual Int_t Filter();
  virtual Int_t DoGenericSubtractionJetMass();
  virtual Int_t DoGenericSubtractionGR(Int_t ijet);
  virtual Int_t DoGenericSubtractionJetAngularity();
  virtual Int_t DoGenericSubtractionJetpTD();
  virtual Int_t DoGenericSubtractionJetCircularity();
  virtual Int_t DoGenericSubtractionJetSigma2();
  virtual Int_t DoGenericSubtractionJetConstituent();
  virtual Int_t DoGenericSubtractionJetLeSub();
  virtual Int_t DoConstituentSubtraction();

  void SetName(const char* name)        { fName           = name;    }
  void SetTitle(const char* title)      { fTitle          = title;   }
  void SetStrategy(const fastjet::Strategy &strat)                 { fStrategy = strat;  }
  void SetAlgorithm(const fastjet::JetAlgorithm &algor)            { fAlgor    = algor;  }
  void SetRecombScheme(const fastjet::RecombinationScheme &scheme) { fScheme   = scheme; }
  void SetAreaType(const fastjet::AreaType &atype)                 { fAreaType = atype;  }
  void SetNRepeats(Int_t nrepeat)       { fNGhostRepeats  = nrepeat; }
  void SetGhostArea(Double_t gharea)    { fGhostArea      = gharea;  }
  void SetMaxRap(Double_t maxrap)       { fMaxRap         = maxrap;  }
  void SetR(Double_t r)                 { fR              = r;       }
  void SetGridScatter(Double_t gridSc)  { fGridScatter    = gridSc;  }
  void SetKtScatter(Double_t ktSc)      { fKtScatter      = ktSc;    }
  void SetMeanGhostKt(Double_t meankt)  { fMeanGhostKt    = meankt;  }
  void SetPluginAlgor(Int_t plugin)     { fPluginAlgor    = plugin;  }
  void SetUseArea4Vector(Bool_t useA4v) { fUseArea4Vector = useA4v;  }
  void SetupAlgorithmfromOpt(const char *option);
  void SetupAreaTypefromOpt(const char *option);
  void SetupSchemefromOpt(const char *option);
  void SetupStrategyfromOpt(const char *option);
  void SetLegacyMode (Bool_t mode)      { fLegacyMode ^= mode; }
  void SetLegacyFJ();
  void SetUseExternalBkg(Bool_t b, Double_t rho, Double_t rhom) { fUseExternalBkg = b; fRho = rho; fRhom = rhom;}
  void SetRMaxAndStep(Double_t rmax, Double_t dr) {fRMax = rmax; fDRStep = dr; }

 protected:
  TString                                fName;               //!
  TString                                fTitle;              //!
  std::vector<fastjet::PseudoJet>        fInputVectors;       //!
  std::vector<fastjet::PseudoJet>        fInputGhosts;        //!
  std::vector<fastjet::PseudoJet>        fInclusiveJets;      //!
  std::vector<fastjet::PseudoJet>        fFilteredJets;       //!
  std::vector<double>                    fSubtractedJetsPt;   //!
  std::vector<fastjet::PseudoJet>        fConstituentSubtrJets; //!
  fastjet::AreaDefinition               *fAreaDef;            //!
  fastjet::VoronoiAreaSpec              *fVorAreaSpec;        //!
  fastjet::GhostedAreaSpec              *fGhostedAreaSpec;    //!
  fastjet::JetDefinition                *fJetDef;             //!
  fastjet::JetDefinition::Plugin        *fPlugin;             //!
  fastjet::Selector                     *fRange;              //!
  fastjet::ClusterSequenceArea          *fClustSeq;           //!
  fastjet::ClusterSequence              *fClustSeqSA;                //!
  fastjet::ClusterSequenceActiveAreaExplicitGhosts *fClustSeqActGhosts; //!
  fastjet::Strategy                      fStrategy;           //!
  fastjet::JetAlgorithm                  fAlgor;              //!
  fastjet::RecombinationScheme           fScheme;             //!
  fastjet::AreaType                      fAreaType;           //!
  Int_t                                  fNGhostRepeats;      //!
  Double_t                               fGhostArea;	      //!
  Double_t                               fMaxRap;	      //!
  Double_t                               fR;                  //!
  // no setters for the moment - used default values in the constructor
  Double_t                               fGridScatter;        //!
  Double_t                               fKtScatter;	      //!
  Double_t                               fMeanGhostKt;        //!
  Int_t                                  fPluginAlgor;        //!
  // extra parameters
  Double_t                               fMedUsedForBgSub;    //!
  Bool_t                                 fUseArea4Vector;     //!

  fastjet::JetMedianBackgroundEstimator   *fBkrdEstimator;    //!
  //from contrib package
  fastjet::contrib::GenericSubtractor     *fGenSubtractor;    //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoJetMass;        //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoGRNum;          //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoGRDen;          //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoJetAngularity;  //!  
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoJetpTD;         //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoJetCircularity; //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoJetSigma2;      //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoJetConstituent; //!
  std::vector<fastjet::contrib::GenericSubtractorInfo> fGenSubtractorInfoJetLeSub;       //!
  
  Bool_t                                   fDoFilterArea;         //!
  Bool_t                                   fLegacyMode;           //!
  Bool_t                                   fUseExternalBkg;       //!
  Double_t                                 fRho;                  //  pT background density
  Double_t                                 fRhom;                 //  mT background density
  Double_t                                 fRMax;             //!
  Double_t                                 fDRStep;           //!
  std::vector<double>                      fGRNumerator;      //!
  std::vector<double>                      fGRDenominator;    //!
  std::vector<double>                      fGRNumeratorSub;   //!
  std::vector<double>                      fGRDenominatorSub; //!

  virtual void   SubtractBackground(const Double_t median_pt = -1);

 private:
  FJWrapper();
  FJWrapper(const FJWrapper& wrapper);
  FJWrapper& operator = (const FJWrapper& wrapper);
};
#endif

#ifdef FJWrapper_cc
#undef FJWrapper_cc 

#if defined __GNUC__
#pragma GCC system_header
#endif

namespace fj = fastjet;

//_________________________________________________________________________________________________
FJWrapper::FJWrapper(const char *name, const char *title)
  : 
    fName              (name)
  , fTitle             (title)
  , fInputVectors      ( )
  , fInputGhosts       ( )
  , fInclusiveJets     ( )
  , fFilteredJets      ( )
  , fSubtractedJetsPt  ( )
  , fConstituentSubtrJets ( )
  , fAreaDef           (0)
  , fVorAreaSpec       (0)
  , fGhostedAreaSpec   (0)
  , fJetDef            (0)
  , fPlugin            (0)
  , fRange             (0)
  , fClustSeq          (0)
  , fClustSeqSA        (0)
  , fClustSeqActGhosts (0)
  , fStrategy          (fj::Best)
  , fAlgor             (fj::kt_algorithm)
  , fScheme            (fj::BIpt_scheme)
  , fAreaType          (fj::active_area)
  , fNGhostRepeats     (1)
  , fGhostArea         (0.005)
  , fMaxRap            (5.)
  , fR                 (0.4)
  , fGridScatter       (1.0)
  , fKtScatter         (0.1)
  , fMeanGhostKt       (1e-100)
  , fPluginAlgor       (0)
  , fMedUsedForBgSub   (0)
  , fUseArea4Vector    (kFALSE)
  , fBkrdEstimator     (0)
  , fGenSubtractor     (0)
  , fGenSubtractorInfoJetMass ( )
  , fGenSubtractorInfoGRNum ( )
  , fGenSubtractorInfoGRDen ( )
  , fGenSubtractorInfoJetAngularity ( )
  , fGenSubtractorInfoJetpTD ( )
  , fGenSubtractorInfoJetCircularity( )
  , fGenSubtractorInfoJetSigma2()
  , fGenSubtractorInfoJetConstituent ( )
  , fGenSubtractorInfoJetLeSub ( )
  , fDoFilterArea      (false)
  , fLegacyMode        (false)
  , fUseExternalBkg    (false)
  , fRho               (0)
  , fRhom              (0)
  , fRMax(2.)
  , fDRStep(0.04)
  , fGRNumerator()
  , fGRDenominator()
  , fGRNumeratorSub()
  , fGRDenominatorSub()
{
  // Constructor.
}

//_________________________________________________________________________________________________
FJWrapper::~FJWrapper()
{  
  // Destructor.

  if (fAreaDef)                 { delete fAreaDef;           fAreaDef         = NULL; }
  if (fVorAreaSpec)             { delete fVorAreaSpec;       fVorAreaSpec     = NULL; }
  if (fGhostedAreaSpec)         { delete fGhostedAreaSpec;   fGhostedAreaSpec = NULL; }
  if (fJetDef)                  { delete fJetDef;            fJetDef          = NULL; }
  if (fPlugin)                  { delete fPlugin;            fPlugin          = NULL; }
  if (fRange)                   { delete fRange;             fRange           = NULL; }
  if (fClustSeq)                { delete fClustSeq;          fClustSeq        = NULL; }
  if (fClustSeqSA)              { delete fClustSeqSA;        fClustSeqSA        = NULL; }
  if (fClustSeqActGhosts)       { delete fClustSeqActGhosts; fClustSeqActGhosts = NULL; }
  if (fBkrdEstimator)           { delete fBkrdEstimator;     fBkrdEstimator   = NULL; }
  if (fGenSubtractor)           { delete fGenSubtractor;     fGenSubtractor   = NULL; }
  // if (fConstituentSubtractor)   { delete fConstituentSubtractor; fConstituentSubtractor = NULL; }

}

//_________________________________________________________________________________________________
void FJWrapper::CopySettingsFrom(const FJWrapper& wrapper)
{
  // Copy some settings.
  // You very often want to keep most of the settings 
  // but change only the algorithm or R - do it after call to this function

  fStrategy         = wrapper.fStrategy;
  fAlgor            = wrapper.fAlgor;
  fScheme           = wrapper.fScheme;
  fAreaType         = wrapper.fAreaType;
  fNGhostRepeats    = wrapper.fNGhostRepeats;
  fGhostArea        = wrapper.fGhostArea;
  fMaxRap           = wrapper.fMaxRap;
  fR                = wrapper.fR;
  fGridScatter      = wrapper.fGridScatter;
  fKtScatter        = wrapper.fKtScatter;
  fMeanGhostKt      = wrapper.fMeanGhostKt;
  fPluginAlgor      = wrapper.fPluginAlgor;
  fUseArea4Vector   = wrapper.fUseArea4Vector;
  fLegacyMode       = wrapper.fLegacyMode;
  fUseExternalBkg   = wrapper.fUseExternalBkg;
  fRho              = wrapper.fRho;
  fRhom             = wrapper.fRhom;
}

//_________________________________________________________________________________________________
void FJWrapper::Clear(const Option_t */*opt*/)
{
  // Simply clear the input vectors.
  // Make sure done on every event if the instance is reused
  // Reset the median to zero.

  fInputVectors.clear();
  fInputGhosts.clear();
  fMedUsedForBgSub = 0;

  // for the moment brute force delete everything
  if (fAreaDef)                 { delete fAreaDef;           fAreaDef         = NULL; }
  if (fVorAreaSpec)             { delete fVorAreaSpec;       fVorAreaSpec     = NULL; }
  if (fGhostedAreaSpec)         { delete fGhostedAreaSpec;   fGhostedAreaSpec = NULL; }
  if (fJetDef)                  { delete fJetDef;            fJetDef          = NULL; }
  if (fPlugin)                  { delete fPlugin;            fPlugin          = NULL; }
  if (fRange)                   { delete fRange;             fRange           = NULL; }
  if (fClustSeq)                { delete fClustSeq;          fClustSeq        = NULL; }
  if (fClustSeqSA)              { delete fClustSeqSA;        fClustSeqSA        = NULL; }
  if (fClustSeqActGhosts)       { delete fClustSeqActGhosts; fClustSeqActGhosts = NULL; }
  if (fBkrdEstimator)           { delete fBkrdEstimator;     fBkrdEstimator   = NULL; }
  if (fGenSubtractor)           { delete fGenSubtractor;     fGenSubtractor   = NULL; }
  // if (fConstituentSubtractor)   { delete fConstituentSubtractor; fConstituentSubtractor = NULL; }
}

//_________________________________________________________________________________________________
void FJWrapper::AddInputVector(Double_t px, Double_t py, Double_t pz, Double_t E, Int_t index)
{
  // Make the input pseudojet.

  fastjet::PseudoJet inVec(px, py, pz, E);
  
  if (index > -99999) {
    inVec.set_user_index(index);
  } else {
    inVec.set_user_index(fInputVectors.size());
  }

  // add to the fj container of input vectors
  fInputVectors.push_back(inVec);
}

//_________________________________________________________________________________________________
void FJWrapper::AddInputVector(const fj::PseudoJet& vec, Int_t index)
{
  // Add an input pseudojet.

  fj::PseudoJet inVec = vec;
  
  if (index > -99999) {
    inVec.set_user_index(index);
  } else {
    inVec.set_user_index(fInputVectors.size());
  }

  // add to the fj container of input vectors
  fInputVectors.push_back(inVec);
}

//_________________________________________________________________________________________________
void FJWrapper::AddInputVectors(const std::vector<fj::PseudoJet>& vecs, Int_t offsetIndex)
{
  // Add the input from vector of pseudojets.

  for (UInt_t i = 0; i < vecs.size(); ++i) {
    fj::PseudoJet inVec = vecs[i];
    if (offsetIndex > -99999)
      inVec.set_user_index(fInputVectors.size() + offsetIndex);
    // add to the fj container of input vectors
    fInputVectors.push_back(inVec);
  }
}

//_________________________________________________________________________________________________
void FJWrapper::AddInputGhost(Double_t px, Double_t py, Double_t pz, Double_t E, Int_t index)
{
  // Make the input pseudojet.

  fastjet::PseudoJet inVec(px, py, pz, E);
  
  if (index > -99999) {
    inVec.set_user_index(index);
  } else {
    inVec.set_user_index(fInputGhosts.size());
  }

  // add to the fj container of input vectors
  fInputGhosts.push_back(inVec);
  if (!fDoFilterArea) fDoFilterArea = kTRUE;
}

//_________________________________________________________________________________________________
Double_t FJWrapper::GetJetArea(UInt_t idx) const
{
  // Get the jet area.

  Double_t retval = -1; // really wrong area..
  if ( idx < fInclusiveJets.size() ) {
    retval = fClustSeq->area(fInclusiveJets[idx]);
  } else {
    Printf("[e] ::GetJetArea wrong index: %d",idx);
  }
  return retval;
}

//_________________________________________________________________________________________________
Double_t FJWrapper::GetFilteredJetArea(UInt_t idx) const
{
  // Get the filtered jet area.

  Double_t retval = -1; // really wrong area..
  if (fDoFilterArea && fClustSeqActGhosts && (idx<fFilteredJets.size())) {
    retval = fClustSeqActGhosts->area(fFilteredJets[idx]);
  } else {
    Printf("[e] ::GetFilteredJetArea wrong index: %d",idx);
  }
  return retval;
}

//_________________________________________________________________________________________________
fastjet::PseudoJet FJWrapper::GetJetAreaVector(UInt_t idx) const
{
  // Get the jet area as vector.
  fastjet::PseudoJet retval;
  if ( idx < fInclusiveJets.size() ) {
    retval = fClustSeq->area_4vector(fInclusiveJets[idx]);
  } else {
    Printf("[e] ::GetJetArea wrong index: %d",idx);
  }
  return retval;
}

//_________________________________________________________________________________________________
fastjet::PseudoJet FJWrapper::GetFilteredJetAreaVector(UInt_t idx) const
{
  // Get the jet area as vector.
  fastjet::PseudoJet retval;
  if (fDoFilterArea && fClustSeqActGhosts && (idx<fFilteredJets.size())) {
    retval = fClustSeqActGhosts->area_4vector(fFilteredJets[idx]);
  } else {
    Printf("[e] ::GetFilteredJetArea wrong index: %d",idx);
  }
  return retval;
}

//_________________________________________________________________________________________________
std::vector<double> FJWrapper::GetSubtractedJetsPts(Double_t median_pt, Bool_t sorted)
{ 
  // Get subtracted jets pTs, returns vector.

  SubtractBackground(median_pt);
  
  if (kTRUE == sorted) {
    std::sort(fSubtractedJetsPt.begin(), fSubtractedJetsPt.begin());
  }
  return fSubtractedJetsPt;
}

//_________________________________________________________________________________________________
Double_t FJWrapper::GetJetSubtractedPt(UInt_t idx) const
{
  // Get subtracted jets pTs, returns Double_t.

  Double_t retval = -99999.; // really wrong pt..
  if ( idx < fSubtractedJetsPt.size() ) {
    retval = fSubtractedJetsPt[idx];
  }
  return retval;
}

//_________________________________________________________________________________________________
std::vector<fastjet::PseudoJet>
FJWrapper::GetJetConstituents(UInt_t idx) const
{
  // Get jets constituents.

  std::vector<fastjet::PseudoJet> retval;
  
  if ( idx < fInclusiveJets.size() ) {
    retval = fClustSeq->constituents(fInclusiveJets[idx]);
  } else {
    Printf("[e] ::GetJetConstituents wrong index: %d",idx);
  }
  
  return retval;
}

//_________________________________________________________________________________________________
std::vector<fastjet::PseudoJet>
FJWrapper::GetFilteredJetConstituents(UInt_t idx) const
{
  // Get jets constituents.

  std::vector<fastjet::PseudoJet> retval;
  
  if ( idx < fFilteredJets.size() ) {
    if (fClustSeqSA)        retval = fClustSeqSA->constituents(fFilteredJets[idx]);
    if (fClustSeqActGhosts) retval = fClustSeqActGhosts->constituents(fFilteredJets[idx]);
  } else {
    Printf("[e] ::GetFilteredJetConstituents wrong index: %d",idx);
  }
  
  return retval;
}

//_________________________________________________________________________________________________
void FJWrapper::GetMedianAndSigma(Double_t &median, Double_t &sigma, Int_t remove) const
{
  // Get the median and sigma from fastjet.
  // User can also do it on his own because the cluster sequence is exposed (via a getter)

  if (!fClustSeq) {
    Printf("[e] Run the jfinder first.");
    return;
  }

  Double_t mean_area = 0;
  try {
    if(0 == remove) {
      fClustSeq->get_median_rho_and_sigma(*fRange, fUseArea4Vector, median, sigma, mean_area);
    }  else {
      std::vector<fastjet::PseudoJet> input_jets = sorted_by_pt(fClustSeq->inclusive_jets());
      input_jets.erase(input_jets.begin(), input_jets.begin() + remove);
      fClustSeq->get_median_rho_and_sigma(input_jets, *fRange, fUseArea4Vector, median, sigma, mean_area);
      input_jets.clear();
    }
  } catch (fj::Error) {
    Printf(" [w] FJ Exception caught.");
    median = -1.;
    sigma = -1;
  }
}

//_________________________________________________________________________________________________
Int_t FJWrapper::Run()
{
  // Run the actual jet finder.

  if (fAreaType == fj::voronoi_area) {
    // Rfact - check dependence - default is 1.
    // NOTE: hardcoded variable!
    fVorAreaSpec = new fj::VoronoiAreaSpec(1.); 
    fAreaDef     = new fj::AreaDefinition(*fVorAreaSpec);      
  } else {
    fGhostedAreaSpec = new fj::GhostedAreaSpec(fMaxRap,
                                               fNGhostRepeats, 
                                               fGhostArea,
                                               fGridScatter,
                                               fKtScatter,
                                               fMeanGhostKt);

    fAreaDef = new fj::AreaDefinition(*fGhostedAreaSpec, fAreaType);
  }
  
  // this is acceptable by fastjet:
  fRange = new fj::Selector(fj::SelectorAbsRapMax(fMaxRap - 0.95 * fR));

  if (fAlgor == fj::plugin_algorithm) {
    if (fPluginAlgor == 0) {
      // SIS CONE ALGOR
      // NOTE: hardcoded split parameter
      Double_t overlap_threshold = 0.75; // NOTE: this actually splits a lot: thr/min(pt1,pt2)
      fPlugin = new fj::SISConePlugin(fR, 
                                      overlap_threshold,
                                      0,    //search of stable cones - zero = until no more
                                      1.0); // this should be seed effectively for proto jets
      fJetDef = new fastjet::JetDefinition(fPlugin);
    } else if (fPluginAlgor == 1) {
      // CDF cone
      // NOTE: hardcoded split parameter
      Double_t overlap_threshold = 0.75; // NOTE: this actually splits a lot: thr/min(pt1,pt2)
      fPlugin = new fj::CDFMidPointPlugin(fR, 
                                      overlap_threshold,
                                      1.0,    //search of stable cones - zero = until no more
                                      1.0); // this should be seed effectively for proto jets
      fJetDef = new fastjet::JetDefinition(fPlugin);
    } else {
      Printf("[e] Unrecognized plugin number!");
    }
  } else {
    fJetDef = new fj::JetDefinition(fAlgor, fR, fScheme, fStrategy);
  }
  
  try {
    fClustSeq = new fj::ClusterSequenceArea(fInputVectors, *fJetDef, *fAreaDef);
  } catch (fj::Error) {
    Printf(" [w] FJ Exception caught.");
    return -1;
  }

  // FJ3 :: Define an JetMedianBackgroundEstimator just in case it will be used 
  fBkrdEstimator     = new fj::JetMedianBackgroundEstimator(fj::SelectorAbsRapMax(fMaxRap));

  if (fLegacyMode) { SetLegacyFJ(); } // for FJ 2.x even if fLegacyMode is set, SetLegacyFJ is dummy

  // inclusive jets:
  fInclusiveJets.clear();
  fInclusiveJets = fClustSeq->inclusive_jets(0.0); 

  return 0;
}

//_________________________________________________________________________________________________
Int_t FJWrapper::Filter()
{
//
//  FJWrapper::Filter
//

  fJetDef = new fj::JetDefinition(fAlgor, fR, fScheme, fStrategy);

  if (fDoFilterArea) {
    if (fInputGhosts.size()>0) {
      try {
        fClustSeqActGhosts = new fj::ClusterSequenceActiveAreaExplicitGhosts(fInputVectors,
                                                                           *fJetDef,
                                                                            fInputGhosts,
                                                                            fGhostArea);
      } catch (fj::Error) {
        Printf(" [w] FJ Exception caught.");
        return -1;
      }

      fFilteredJets.clear();
      fFilteredJets =  fClustSeqActGhosts->inclusive_jets(0.0); 
    } else {
      return -1;
    }
  } else {
    try {
      fClustSeqSA = new fastjet::ClusterSequence(fInputVectors, *fJetDef);
    } catch (fj::Error) {
      Printf(" [w] FJ Exception caught.");
      return -1;
    }

    fFilteredJets.clear();
    fFilteredJets = fClustSeqSA->inclusive_jets(0.0);
  }

  return 0;
}

//_________________________________________________________________________________________________
void FJWrapper::SetLegacyFJ()
{
  // This methods enable legacy behaviour (FastJet 2.x) when AliROOT is compiled with FastJet 3.x
  std::cout << "WARNING! Setting FastJet in legacy mode" << std::endl;
  if (fGhostedAreaSpec) { fGhostedAreaSpec->set_fj2_placement(kTRUE); }
  if (fBkrdEstimator) {
    fBkrdEstimator->set_provide_fj2_sigma(kTRUE);
    fBkrdEstimator->set_use_area_4vector(kFALSE);
  } 
}

//_________________________________________________________________________________________________
void FJWrapper::SubtractBackground(Double_t median_pt)
{
  // Subtract the background (specify the value - see below the meaning).
  // Negative argument means the bg will be determined with the current algorithm
  // this is the default behavior. Zero is allowed
  // Note: user may set the switch for area4vector based subtraction.

  Double_t median    = 0;
  Double_t sigma     = 0;
  Double_t mean_area = 0;

  // clear the subtracted jet pt's vector<double>
  fSubtractedJetsPt.clear();

  // check what was specified (default is -1)
  if (median_pt < 0) {
    try {
      fClustSeq->get_median_rho_and_sigma(*fRange, fUseArea4Vector, median, sigma, mean_area);
    }

    catch (fj::Error) {
      Printf(" [w] FJ Exception caught.");
      median = -9999.;
      sigma = -1;
      fMedUsedForBgSub = median;
      return;
    }
    fMedUsedForBgSub = median;
  } else {
    // we do not know the sigma in this case
    sigma = -1;
    if (0.0 == median_pt) {
      fMedUsedForBgSub = 0.;
    } else {
      fMedUsedForBgSub = median_pt;
    }
  }

  // subtract:
  for (unsigned i = 0; i < fInclusiveJets.size(); i++) {
    if ( fUseArea4Vector ) {
      // subtract the background using the area4vector
      fj::PseudoJet area4v = fClustSeq->area_4vector(fInclusiveJets[i]);
      fj::PseudoJet jet_sub = fInclusiveJets[i] - area4v * fMedUsedForBgSub;
      fSubtractedJetsPt.push_back(jet_sub.perp()); // here we put only the pt of the jet - note: this can be negative
    } else {
      // subtract the background using scalars
      // fj::PseudoJet jet_sub = fInclusiveJets[i] - area * fMedUsedForBgSub_;
      Double_t area = fClustSeq->area(fInclusiveJets[i]);
      // standard subtraction
      Double_t pt_sub = fInclusiveJets[i].perp() - fMedUsedForBgSub * area;
      fSubtractedJetsPt.push_back(pt_sub); // here we put only the pt of the jet - note: this can be negative
    }
  }
}

//_________________________________________________________________________________________________
Int_t FJWrapper::DoGenericSubtractionJetMass() {
  
  //Do generic subtraction for jet mass
  if(fUseExternalBkg)   fGenSubtractor     = new fj::contrib::GenericSubtractor(fRho,fRhom);
  else                  fGenSubtractor     = new fj::contrib::GenericSubtractor(fBkrdEstimator);

  // Define jet shape
  FJJetShapeMass shapeMass;

  // clear the generic subtractor info vector
  fGenSubtractorInfoJetMass.clear();
  for (unsigned i = 0; i < fInclusiveJets.size(); i++) {
    fj::contrib::GenericSubtractorInfo info;
    if(fInclusiveJets[i].perp()>1.e-4)
      double subtracted_shape = (*fGenSubtractor)(shapeMass, fInclusiveJets[i], info);
    fGenSubtractorInfoJetMass.push_back(info);
  }
  
  return 0;
}

//_________________________________________________________________________________________________
Int_t FJWrapper::DoGenericSubtractionGR(Int_t ijet) {
  
  //Do generic subtraction for jet structure function
  if(fUseExternalBkg)   fGenSubtractor     = new fj::contrib::GenericSubtractor(fRho,fRhom);
  else                  fGenSubtractor     = new fj::contrib::GenericSubtractor(fBkrdEstimator);

  if(ijet>fInclusiveJets.size()) return 0;

  fGRNumerator.clear();
  fGRDenominator.clear();
  fGRNumeratorSub.clear();
  fGRDenominatorSub.clear();

  // Define jet shape
  for(Double_t r = 0.; r<fRMax; r+=fDRStep) {
    FJJetShapeGRNum shapeGRNum(r,fDRStep);
    FJJetShapeGRDen shapeGRDen(r,fDRStep);

    // clear the generic subtractor info vector
    fGenSubtractorInfoGRNum.clear();
    fGenSubtractorInfoGRDen.clear();
    fj::contrib::GenericSubtractorInfo infoNum;
    fj::contrib::GenericSubtractorInfo infoDen;
    if(fInclusiveJets[ijet].perp()>1.e-4) {
      double sub_num = (*fGenSubtractor)(shapeGRNum, fInclusiveJets[ijet], infoNum);
      double sub_den = (*fGenSubtractor)(shapeGRDen, fInclusiveJets[ijet], infoDen);
    }
    fGenSubtractorInfoGRNum.push_back(infoNum);
    fGenSubtractorInfoGRDen.push_back(infoDen);
    fGRNumerator.push_back(infoNum.unsubtracted());
    fGRDenominator.push_back(infoDen.unsubtracted());
    fGRNumeratorSub.push_back(infoNum.second_order_subtracted());
    fGRDenominatorSub.push_back(infoDen.second_order_subtracted());
  }
  
  return 0;
}
//_________________________________________________________________________________________________
Int_t FJWrapper::DoGenericSubtractionJetAngularity() {
  
  //Do generic subtraction for jet mass
  if(fUseExternalBkg)   fGenSubtractor     = new fj::contrib::GenericSubtractor(fRho,fRhom);
  else                  fGenSubtractor     = new fj::contrib::GenericSubtractor(fBkrdEstimator);

  // Define jet shape
  FJJetShapeAngularity shapeAngularity;

  // clear the generic subtractor info vector
  fGenSubtractorInfoJetAngularity.clear();
  for (unsigned i = 0; i < fInclusiveJets.size(); i++) {
    fj::contrib::GenericSubtractorInfo infoAng;
    if(fInclusiveJets[i].perp()>1.e-4)
      double subtracted_shape = (*fGenSubtractor)(shapeAngularity, fInclusiveJets[i], infoAng);
    fGenSubtractorInfoJetAngularity.push_back(infoAng);
  }
  
  return 0;
}
//_________________________________________________________________________________________________
Int_t FJWrapper::DoGenericSubtractionJetpTD() {
  
  //Do generic subtraction for pTD
  if(fUseExternalBkg)   fGenSubtractor     = new fj::contrib::GenericSubtractor(fRho,fRhom);
  else                  fGenSubtractor     = new fj::contrib::GenericSubtractor(fBkrdEstimator);

  // Define jet shape
  FJJetShapepTD shapepTD;

  // clear the generic subtractor info vector
  fGenSubtractorInfoJetpTD.clear();
  for (unsigned i = 0; i < fInclusiveJets.size(); i++) {
    fj::contrib::GenericSubtractorInfo infopTD;
    if(fInclusiveJets[i].perp()>1.e-4)
      double subtracted_shape = (*fGenSubtractor)(shapepTD, fInclusiveJets[i], infopTD);
    fGenSubtractorInfoJetpTD.push_back(infopTD);
  }
  
  return 0;
}
//_________________________________________________________________________________________________
Int_t FJWrapper::DoGenericSubtractionJetCircularity() {
  
  //Do generic subtraction for jet circularity
  if(fUseExternalBkg)   fGenSubtractor     = new fj::contrib::GenericSubtractor(fRho,fRhom);
  else                  fGenSubtractor     = new fj::contrib::GenericSubtractor(fBkrdEstimator);

  // Define jet shape
  FJJetShapeCircularity shapecircularity;

  // clear the generic subtractor info vector
  fGenSubtractorInfoJetCircularity.clear();
  for (unsigned i = 0; i < fInclusiveJets.size(); i++) {
    fj::contrib::GenericSubtractorInfo infoCirc;
    if(fInclusiveJets[i].perp()>1.e-4)
      double subtracted_shape = (*fGenSubtractor)(shapecircularity, fInclusiveJets[i], infoCirc);
    fGenSubtractorInfoJetCircularity.push_back(infoCirc);
  }
  
 return 0;
}
//_________________________________________________________________________________________________
Int_t FJWrapper::DoGenericSubtractionJetSigma2() {
  
  //Do generic subtraction for jet sigma2
  if(fUseExternalBkg)   fGenSubtractor     = new fj::contrib::GenericSubtractor(fRho,fRhom);
  else                  fGenSubtractor     = new fj::contrib::GenericSubtractor(fBkrdEstimator);

  // Define jet shape
  FJJetShapeSigma2 shapesigma2;

  // clear the generic subtractor info vector
  fGenSubtractorInfoJetSigma2.clear();
  for (unsigned i = 0; i < fInclusiveJets.size(); i++) {
    fj::contrib::GenericSubtractorInfo infoSigma;
    if(fInclusiveJets[i].perp()>1.e-4)
      double subtracted_shape = (*fGenSubtractor)(shapesigma2, fInclusiveJets[i], infoSigma);
    fGenSubtractorInfoJetSigma2.push_back(infoSigma);
  }
  
  return 0;
}
//_________________________________________________________________________________________________
Int_t FJWrapper::DoGenericSubtractionJetConstituent() {
  
  //Do generic subtraction for #constituents
  if(fUseExternalBkg)   fGenSubtractor     = new fj::contrib::GenericSubtractor(fRho,fRhom);
  else                  fGenSubtractor     = new fj::contrib::GenericSubtractor(fBkrdEstimator);

  // Define jet shape
  FJJetShapeConstituent shapeconst;

  // clear the generic subtractor info vector
   fGenSubtractorInfoJetConstituent.clear();
  for (unsigned i = 0; i < fInclusiveJets.size(); i++) {
    fj::contrib::GenericSubtractorInfo infoConst;
    if(fInclusiveJets[i].perp()>1.e-4)
      double subtracted_shape = (*fGenSubtractor)(shapeconst, fInclusiveJets[i], infoConst);
    fGenSubtractorInfoJetConstituent.push_back(infoConst);
  }
  
  return 0;
}

//_________________________________________________________________________________________________
Int_t FJWrapper::DoGenericSubtractionJetLeSub() {
  
  //Do generic subtraction for leading-subleading constituent
  if(fUseExternalBkg)   fGenSubtractor     = new fj::contrib::GenericSubtractor(fRho,fRhom);
  else                  fGenSubtractor     = new fj::contrib::GenericSubtractor(fBkrdEstimator);

  // Define jet shape
  FJJetShapeLeSub shapeLeSub;

  // clear the generic subtractor info vector
  fGenSubtractorInfoJetLeSub.clear();
  for (unsigned i = 0; i < fInclusiveJets.size(); i++) {
    fj::contrib::GenericSubtractorInfo infoLeSub;
    if(fInclusiveJets[i].perp()>1.e-4)
      double subtracted_shape = (*fGenSubtractor)(shapeLeSub, fInclusiveJets[i], infoLeSub);
    fGenSubtractorInfoJetLeSub.push_back(infoLeSub);
  }
  
  return 0;
}

//_________________________________________________________________________________________________
Int_t FJWrapper::DoConstituentSubtraction() {
  //Do constituent subtraction
  fj::contrib::ConstituentSubtractor *subtractor;
  if(fUseExternalBkg)
    subtractor     = new fj::contrib::ConstituentSubtractor(fRho,fRhom,0.,-1.);//kFALSE,kTRUE);
  else                 subtractor     = new fj::contrib::ConstituentSubtractor(fBkrdEstimator);

  //clear constituent subtracted jets
  fConstituentSubtrJets.clear();
  for (unsigned i = 0; i < fInclusiveJets.size(); i++) {
    fj::PseudoJet subtracted_jet(0.,0.,0.,0.);
    if(fInclusiveJets[i].perp()>0.)
      subtracted_jet = (*subtractor)(fInclusiveJets[i]);
    fConstituentSubtrJets.push_back(subtracted_jet);
  }
  if(subtractor) {delete subtractor; subtractor = 0; }
  
  return 0;
}

//_________________________________________________________________________________________________
void FJWrapper::SetupAlgorithmfromOpt(const char *option)
{
  // Setup algorithm from char.

  std::string opt(option);
  
  if (!opt.compare("kt"))                fAlgor    = fj::kt_algorithm;
  if (!opt.compare("antikt"))            fAlgor    = fj::antikt_algorithm;
  if (!opt.compare("cambridge"))         fAlgor    = fj::cambridge_algorithm;
  if (!opt.compare("genkt"))             fAlgor    = fj::genkt_algorithm;
  if (!opt.compare("cambridge_passive")) fAlgor    = fj::cambridge_for_passive_algorithm;
  if (!opt.compare("genkt_passive"))     fAlgor    = fj::genkt_for_passive_algorithm;
  if (!opt.compare("ee_kt"))             fAlgor    = fj::ee_kt_algorithm;
  if (!opt.compare("ee_genkt"))          fAlgor    = fj::ee_genkt_algorithm;
  if (!opt.compare("plugin"))            fAlgor    = fj::plugin_algorithm;
}

//_________________________________________________________________________________________________
void FJWrapper::SetupAreaTypefromOpt(const char *option)
{
  // Setup area type from char.

  std::string opt(option);

  if (!opt.compare("active"))                      fAreaType = fj::active_area;
  if (!opt.compare("invalid"))                     fAreaType = fj::invalid_area;
  if (!opt.compare("active_area_explicit_ghosts")) fAreaType = fj::active_area_explicit_ghosts;
  if (!opt.compare("one_ghost_passive"))           fAreaType = fj::one_ghost_passive_area;
  if (!opt.compare("passive"))                     fAreaType = fj::passive_area;
  if (!opt.compare("voronoi"))                     fAreaType = fj::voronoi_area;
}

//_________________________________________________________________________________________________
void FJWrapper::SetupSchemefromOpt(const char *option)
{
  //
  // setup scheme from char
  //

  std::string opt(option);

  if (!opt.compare("BIpt"))   fScheme   = fj::BIpt_scheme;
  if (!opt.compare("BIpt2"))  fScheme   = fj::BIpt2_scheme;
  if (!opt.compare("E"))      fScheme   = fj::E_scheme;
  if (!opt.compare("pt"))     fScheme   = fj::pt_scheme;
  if (!opt.compare("pt2"))    fScheme   = fj::pt2_scheme;
  if (!opt.compare("Et"))     fScheme   = fj::Et_scheme;
  if (!opt.compare("Et2"))    fScheme   = fj::Et2_scheme;
}

//_________________________________________________________________________________________________
void FJWrapper::SetupStrategyfromOpt(const char *option)
{
  // Setup strategy from char.

  std::string opt(option);
  
  if (!opt.compare("Best"))            fStrategy = fj::Best;
  if (!opt.compare("N2MinHeapTiled"))  fStrategy = fj::N2MinHeapTiled;
  if (!opt.compare("N2Tiled"))         fStrategy = fj::N2Tiled;
  if (!opt.compare("N2PoorTiled"))     fStrategy = fj::N2PoorTiled;
  if (!opt.compare("N2Plain"))         fStrategy = fj::N2Plain;
  if (!opt.compare("N3Dumb"))          fStrategy = fj::N3Dumb;
  if (!opt.compare("NlnN"))            fStrategy = fj::NlnN;
  if (!opt.compare("NlnN3pi"))         fStrategy = fj::NlnN3pi;
  if (!opt.compare("NlnN4pi"))         fStrategy = fj::NlnN4pi;
  if (!opt.compare("NlnNCam4pi"))      fStrategy = fj::NlnNCam4pi;
  if (!opt.compare("NlnNCam2pi2R"))    fStrategy = fj::NlnNCam2pi2R;
  if (!opt.compare("NlnNCam"))         fStrategy = fj::NlnNCam;
  if (!opt.compare("plugin"))          fStrategy = fj::plugin_strategy;
}
#endif
